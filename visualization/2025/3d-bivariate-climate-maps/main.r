# 1. PACKAGES
#------------------

install.packages("remotes")
remotes::install_github(
    "inSileco/rchelsa"
)
remotes::install_github(
    "chris-prener/biscale"
)

# install and load all packages
pacman::p_load(
    geodata, tidyverse, sf, terra,
    rchelsa, biscale, elevatr, cowplot,
    gridGraphics, rayshader
)

# 2. CHELSA DATA
#----------------

# set the working directory
main_dir <- getwd()

# define a vector of IDs to download
ids <- c(1, 12)

# function to download CHELSA data
download_chelsa_data <- function(id, path){
    rchelsa::get_chelsea_data(
        categ = "clim", type = "bio",
        id = id, path = path
    )
}

# download data for each id
lapply(ids, download_chelsa_data, path = main_dir)

list.files()

# load the raster files
temp <- terra::rast("CHELSA_bio10_01.tif")
prec <- terra::rast("CHELSA_bio10_12.tif")

# average precipitation
prec_average <- prec / 30

# Combine average temperature and precipitation
# into a raster stack
temp_prec <- c(temp, prec_average)

# assign names to each layer in the stack
names(temp_prec) <- c("temperature", "precipitation")

# 3. COUNTRY POLYGON
#-------------------

country_sf <- geodata::gadm(
    country = "BEL", level = 0,
    path = main_dir
) |>
sf::st_as_sf()

# 4. CROP AND RESAMPLE
#---------------------

# define the target CRS
target_crs <- "EPSG:3035"

# crop the input raster to the
# country's extent and apply a mask
temp_prec_country <- terra::crop(
    temp_prec, country_sf,
    mask = TRUE
)

# Obtain AWS tiles DEM data from elevatr
# convert to terra SpatRaster and crop
dem <- elevatr::get_elev_raster(
    locations = country_sf, z = 8,
    clip = "locations"
) |> terra::rast() |>
terra::crop(country_sf, mask = TRUE)

# resample the raster to match DEM resolution
# using bilinear interpolation, then reproject

temp_prec_resampled <- terra::resample(
    x = temp_prec_country,
    y = dem, method = "bilinear"
) |> terra::project(target_crs)

# plot the resampled raster
terra::plot(temp_prec_resampled)

# convert the raster to dataframe with coordinates
temp_prec_df <- as.data.frame(
    temp_prec_resampled, xy = TRUE
)

# 5. BREAKS, PALETTE AND PLOT THEME
#----------------------------------

# create bivariate classes using biscale
breaks <- biscale::bi_class(
    temp_prec_df, x = temperature,
    y = precipitation, style = "fisher",
    dim = 3
)

# Define the color palette
pal <- "DkBlue"

# define a custom theme for the map
theme_for_the_win <- function(){
    theme_minimal() +
    theme(
        axis.title = element_blank(),
        plot.background = element_rect(
            fill = "white", color = NA
        ),
        plot.title = element_text(
            color = "grey10", hjust = .5,
            face = "bold", vjust = -1
        ),
        plot.subtitle = element_text(
            hjust = .5, vjust = -1
        ),
        plot.caption = element_text(
            size = 9, color = "grey20",
            hjust = .5, vjust = 1
        ),
        plot.margin = unit(c(0, 0, 0, 0), "lines"
        )
    )
}

# 6. 2D BIVARIATE MAP
#--------------------

# create the bivariate map using ggplot2
map <- ggplot(breaks) +
    geom_raster(
        aes(
            x = x, y = y, fill = bi_class
        ), show.legend = TRUE # FALSE
    ) +
    biscale::bi_scale_fill(
        pal = pal, dim = 3,
        flip_axes = TRUE, rotate_pal = FALSE
    ) +
    labs(
        title = "BELGIUM: Temperature and Precipitation",
        subtitle = "Average temperature and precipitation (1981-2010)",
        caption = "Source: CHELSA | Author: Milos makes maps",
        x = "", y = ""
    ) +
    coord_sf(crs = target_crs) +
    theme_for_the_win()

# create the legend for the bivariate map
legend <- biscale::bi_legend(
    pal = pal,
    flip_axes = TRUE,
    rotate_pal = FALSE,
    dim = 3,
    xlab = "Temperature (°C)",
    ylab = "Precipitation (mm)",
    size = 8
)

# combine the map and legend using cowplot
full_map <- cowplot::ggdraw() +
    cowplot::draw_plot(
        plot = map, x = 0, y = 0,
        width = 1, height = 1
    ) +
    cowplot::draw_plot(
        plot = legend, x = .05, y = .13,
        width = .25, height = .25
    )

# display the final map with legend
print(full_map)

# save as PNG file
ggsave(
    filename = "belgium_bivariate_2d.png",
    width = 7, height = 7, dpi = 600,
    device = "png", bg = "white", full_map
)

# 7. CREATE TERRAIN LAYER
#------------------------

# project and convert to DEM to dataframe
dem_df <- dem |>
    terra::project(target_crs) |>
    as.data.frame(xy = TRUE, na.rm = TRUE)

# rename the third column to "dem"
names(dem_df)[3] <- "dem"

# create the terrain layer map
dem_map <- ggplot(
    dem_df, aes(x = x, y = y, fill = dem)
) +
geom_raster() +
scale_fill_gradientn(colors = "white") +
guides(fill = "none") +
labs(
    title = "BELGIUM: Temperature and Precipitation",
    subtitle = "Average temperature and precipitation (1981-2010)",
    caption = "Source: CHELSA | Author: Milos makes maps"
) +
coord_sf(crs = target_crs) +
theme_for_the_win() +
theme(legend.position = "none")

# 8. RENDER 3D SCENE
#-------------------

rayshader::plot_gg(
    ggobj = full_map,
    ggobj_height = dem_map,
    width = 7,
    height = 7,
    windowsize = c(600, 600),
    scale = 100,
    shadow = TRUE,
    shadow_intensity = 1,
    phi = 87, theta = 0, zoom = .56,
    multicore = TRUE
)

# zoom out
rayshader::render_camera(zoom = .6)

# 9. LIGHTS
#----------

url <- "https://dl.polyhaven.org/file/ph-assets/HDRIs/hdr/4k/brown_photostudio_02_4k.hdr"
hdri_file <- basename(url)

download.file(
    url = url,
    destfile = hdri_file,
    mode = "wb"
)

# 10. RENDER 3D OBJECT
#---------------------

rayshader::render_highquality(
    filename = "belgium-bivariate-3d.png",
    preview = TRUE,
    light = FALSE,
    environment_light = hdri_file,
    intensity = 1,
    rotate_env = 90,
    parallel = TRUE,
    width = 2000, height = 2000,
    interactive = FALSE
)
