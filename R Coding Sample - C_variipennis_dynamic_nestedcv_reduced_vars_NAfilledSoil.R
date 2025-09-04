#' ---
#' title: "Culicoides SDM Project: Maxent Modeling Workflow Using SDMtune Package"
#' author: "Melanie Veron"
#' ---
###
###
###
## Culicoides SDM Project: Maxent Modeling Workflow Using SDMtune Package
##
## Species: Culicoides variipennis
##
## Environmental Data Type and Source:
## - Topology: EarthEnv - Elevation (used only as reference layer)
## - Climate: ClimateNA - Annual and seasonal decadal normals (1911-2020) and 30-year normal (1991-2020)
## - Hydrology: HydroSHEDS - Proximity to rivers and lakes
## - Soil: SoilGrids - Bulk density, soil composition (clay/sand/silt), coarse fragments, organic carbon content, pH in water, nitrogen, and volumetric water content (10, 33, and 1500 kPa)
## - Land Cover % Consensus: EarthEnv - 10 land cover classes
##
## * Environmental Raster Harmonization Baseline Parameters:
## - Geographic CRS (WSG84 datum) with latitude-longitude coordinates in degrees
## - Pixel/cell resolution of 0.050908 degrees squared (approximately equivalent to 4 km squared based on 45 degrees latitude)
## - Map extent from southern Canada through Mexico (equal to extent of study areas combined from Culicoides and Endemic Mexico VSV occurrences)
##
## Workflow Parameters:
##
## * Study Area:
## - Calibration and Baseline Prediction: Polygon mask of intersecting and adjacent ecoregions (Ecoregions2017) relative to the occurrence data (subject to availability of environmental data within the clipped area; note that static version of environmental variables used for prediction)
## - Expanded Prediction: Bounding box of the calibration area (subject to availability of environmental data within the cropped area; note that static version of environmental variables used for prediction)
##
## * Excludes occurrence data from Schmidtmann et al. 2011
##
## * Nested cross-validation (to obtain estimated generalization error of modeling procedure):
## - Outer loop (evaluation): Random 5-fold (k = 5) cross-validation
## - Inner loop (calibration): Checkerboard 4-fold (k = 4) cross-validation
##
## * Includes reduced list of variables shared by all Culicoides study species, using variables with VIF < 10 or combination of VIF < 22 & jackknife AUC > 0.6
##
## * Dynamic SDM
###
###
###


# Prepare data for modeling and analysis ----

### Set working environment ----

# prevent "java.lang.OutOfMemoryError: Java heap space" error when running Maxent
if (Sys.info()[["sysname"]] == "Windows") {
  options(java.parameters = "-Xmx1g")
  gc()
} else {
  options(java.parameters = "-Xmx30g")
  gc()
}

script_start_time <- proc.time()

# packages used in this script:
# - tidyverse: to clean, manipulate, and visualize data
# - dismo: to handle various species distribution modeling tasks
# - spThin: to spatially filter occurrence data
# - SDMtune: to train and evaluate species distribution models
# - ENMeval: to partition data into cross-validation folds
# - usdm: to deal with multicollinearity
# - ecospat: to evaluate species distribution models
# - zeallot: to unpack training/testing data set assignment
# - raster: to handle spatial data
# - terra: to handle spatial data
# - sf: to handle spatial data
# - tmap: to plot spatial data
# - ragg: to save plots
# - rasterVis: to plot raster objects
# - viridis: to customize plot color palettes
# - pryr: to save lines of code as pseudo-objects (useful for saving plots)
# - doParallel: to run code in parallel
# - doFuture: to run code in parallel
# - future.apply: to run code in parallel

# install (if necessary) and load packages

if ("pacman" %in% rownames(installed.packages()) == FALSE) {
  install.packages("pacman")
}
library("pacman")
library("devtools")

if (Sys.info()[["sysname"]] == "Windows") {
  options(devtools.path = "C:/Users/melanie/AppData/Local/R/R-dev")
} else {
  options(devtools.path = "/home/melanie.veron/R/R-dev")
}

if (p_isinstalled(SDMtune) == FALSE |
    (p_isinstalled(SDMtune) == TRUE &&
     p_version(SDMtune) != "1.2.0")) {
  devtools::install_version("SDMtune", "1.2.0")
}

{
  devtools::dev_mode()
  {
    p_unload(SDMtune)
    p_unlock()
    if (p_isinstalled(SDMtune) == FALSE |
        (p_isinstalled(SDMtune) == TRUE &&
         p_version(SDMtune) != "1.1.6")) {
      devtools::install_version("SDMtune", "1.1.6")
    }
    p_unlock()
    library(SDMtune, lib.loc = .libPaths()[1])
  }
} # turn dev_mode on and install/load older version (1.1.6) of SDMtune package to prevent errors/crashing when using certain SDMtune functions; as needed, run the problematic SDMtune code with these dev_mode code chunks framing it

{
  p_unload(SDMtune)
  devtools::dev_mode()
  p_load(SDMtune)
} # turn dev_mode off and revert back to SDMtune package version 1.2.0

p_load(
  tidyverse,
  spThin,
  dismo,
  SDMtune,
  ENMeval,
  usdm,
  ecospat,
  plotrix,
  zeallot,
  raster,
  terra,
  sf,
  tmap,
  ragg,
  rasterVis,
  viridis,
  pryr,
  doParallel,
  doFuture,
  future.apply
)

# set up parallel coding
ncores <- availableCores() / 2
if (Sys.info()[["sysname"]] == "Windows" |
    Sys.info()[["sysname"]] == "Mac") {
  num_parallel_workers <- ncores - 3
} else {
  num_parallel_workers <- ncores
}

# get citation and version information about loaded packages
sink(
  "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/loaded_packages_info.txt"
)
print(
  loaded_packages_info <- list(
    citation = sapply(p_loaded(), p_cite, copy2clip = FALSE),
    version = sapply(p_loaded(), p_ver, simplify = FALSE)
  ) %>% list_transpose()
)
sink()

# define absolute working directory (that includes environmental data directory)
if (Sys.info()[["sysname"]] == "Windows") {
  main_path <- "D:/Documents/WORK/Maxent/"
} else {
  main_path <-
    "/project/disease_ecology/melanie.veron/vs_vectors_sdm/"
}

# create and assign temp directory for raster objects
if (dir.exists(
  "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/raster_tempdir"
) == FALSE) {
  dir.create(
    "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/raster_tempdir"
  )
}
raster_tempdir <-
  "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/raster_tempdir/"

# define custom color palette(s)

# source: https://coolors.co/palette/f94144-f3722c-f8961e-f9c74f-90be6d-43aa8b-577590
custom_color_pal_1 <-
  tibble(
    red_imperial = "#f94144",
    orange1_crayola = "#f3722c",
    orange2_carrot = "#f8961e",
    yellow_saffron = "#f9c74f",
    green1_pistachio = "#90be6d",
    green2_zomp = "#43aa8b",
    gray_paynes = "#577590"
  )


### Prepare presence locations ----

# import and check species occurrence data
occurrence <-
  read.csv("data/occurrence_data/occurrence_new_3.csv") %>% dplyr::select(1:22)
head(occurrence)
str(occurrence)
dim(occurrence)

# create presence data frame filtered by species
presence_data <-
  occurrence %>% filter(species == "Culicoides variipennis") %>%
  filter(source != "Data from Schmidtmann et al. 2011")
head(presence_data)
dim(presence_data)

# save pre-thinned presence location data for later comparison
presence_unthinned <-
  presence_data %>% dplyr::select(longitude, latitude)

# perform spatial thinning for the species occurrence records
unlink(
  "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/presence_thinned_list_log.txt"
)
set.seed(25)
presence_thinned_list <- thin(
  loc.data = presence_data,
  lat.col = "latitude",
  long.col = "longitude",
  spec.col = "species",
  thin.par = 10,
  reps = 100,
  locs.thinned.list.return = TRUE,
  write.files = FALSE,
  write.log.file = TRUE,
  log.file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/presence_thinned_list_log.txt"
)

# order list of randomly thinned presence data frames so that the ones with the most records are listed first (in case all the data frames in the list don't have the same number of rows)
presence_thinned_list <-
  presence_thinned_list[order(sapply(presence_thinned_list, nrow), decreasing = TRUE)]

# select one of the thinned data sets that has the max number of records after thinning (per the output) and use it to overwrite the presence data frame; our goal is to filter the input species presence data by the spatial thinning results in order to retain all the information/columns for each record that made it through the thinning process
presence_thinned <- presence_thinned_list[[1]]
dim(presence_thinned)
presence_thinned_coords <-
  presence_thinned %>% as.data.frame %>% rename(longitude = Longitude, latitude = Latitude)

# turn the row names into a column for each data frame to be used as a key for a table join
presence_data <-
  rownames_to_column(presence_data, "row_names")
presence_thinned_coords <-
  rownames_to_column(presence_thinned_coords, "row_names")

# create a new column grouping the record collection year by decade
presence_data <-
  presence_data %>% mutate(decade = floor(year / 10) * 10) %>% drop_na(decade) %>% filter(decade >= 1910)  # make sure to exclude any records collected before decade of 1910

# filter the full presence data frame so that only the thinned records remain
presence_data_filtered <-
  semi_join(presence_data,
            presence_thinned_coords,
            by = "row_names")
head(presence_data_filtered)
dim(presence_data_filtered)

# store the pre-thinned presence data into a new object
presence_data_unfiltered <- presence_data

# plot a bar graph showing the number of species presence records per collection decade (before and after spatial filtering)
presence_data_unfiltered_and_filtered <-
  bind_rows(
    list("Unfiltered" = presence_data_unfiltered, "Filtered" = presence_data_filtered),
    .id = "type"
  )

presence_data_unfiltered_and_filtered$type <-
  factor(presence_data_unfiltered_and_filtered$type,
         levels = presence_data_unfiltered_and_filtered$type %>% unique())

(
  presence_data_unfiltered_and_filtered_timestep_bargraph <-
    ggplot(
      presence_data_unfiltered_and_filtered,
      aes(x = factor(decade, levels = seq(1910, 2010, 10)),
          fill = type)
    ) +
    geom_bar(position = "dodge", color = "black") +
    labs(x = "Decade",
         y = "Number of presence records collected",
         fill = "") +
    scale_x_discrete(drop = FALSE) +
    geom_text(
      stat = "count",
      aes(label = after_stat(count)),
      vjust = -0.4,
      position = position_dodge(0.9)
    ) +
    scale_fill_viridis(
      option = "mako",
      begin = 0.6,
      end = 0.9,
      discrete = TRUE
    ) +
    ggpubr::theme_pubr() +
    theme(axis.title.x = element_text(
      size = 14,
      margin = margin(10, 10, 10, 10)
    )) +
    theme(axis.title.y = element_text(
      size = 14,
      margin = margin(10, 10, 10, 10)
    )) +
    theme(axis.text = element_text(size = 10))
)

# save the image output
pdf(
  "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/presence_data_unfiltered_and_filtered_timestep_bargraph.pdf",
  width = 10,
  height = 6
)
presence_data_unfiltered_and_filtered_timestep_bargraph
dev.off()
agg_png(
  "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/presence_data_unfiltered_and_filtered_timestep_bargraph.png",
  width = 10,
  height = 6,
  units = "in",
  res = 400
)
presence_data_unfiltered_and_filtered_timestep_bargraph
dev.off()

# get the frequency distribution table for species presence records by collection decade
(freq_table_timestep <- table(presence_data_filtered$decade))

# save output as a csv file
write.csv(freq_table_timestep,
          file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/freq_table_timestep.csv",
          row.names = FALSE)

# get the proportion distribution table for species presence records by collection decade
num_observations <- nrow(presence_data_filtered)
(prob_table_timestep <-
    freq_table_timestep / num_observations)

# save output as a csv file
write.csv(prob_table_timestep,
          file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/prob_table_timestep.csv",
          row.names = FALSE)

# get just the coordinate values from the filtered presence data set
presence <-
  presence_data_filtered %>% dplyr::select(longitude, latitude)
head(presence)


### Acquire environmental variables ----

# import raster files (environmental variables/predictors)
files <-
  list.files(
    c(
      str_c(
        main_path,
        "environmental_data/NEW_environmental_data/ClimateNA/final/normal_1991_2020/Normal_1991_2020MSY"
      ),
      str_c(
        main_path,
        "environmental_data/NEW_environmental_data/EarthEnv_Consensus_Land_Cover/final"
      ),
      str_c(
        main_path,
        "environmental_data/NEW_environmental_data/HydroSHEDS_Hydrology/final"
      ),
      str_c(
        main_path,
        "environmental_data/NEW_environmental_data/SoilGrids/final_NAfilled"
      )
    ),
    pattern = 'tif$',
    full.names = TRUE,
    recursive = FALSE
  )
files <- setdiff(files, str_subset(files, "monthly"))
head(files)
tail(files)
length(files)

# convert the imported files into a rasterstack object
env <- raster::stack(files) %>% raster::trim()

# check the env rasterstack object
env

# check the environmental variable names
names(env)

# import spreadsheet containing readable variable names, units, and other metadata -- then filter by env layer names specific to this script
env_layer_names <-
  read.csv(
    str_c(
      main_path,
      "environmental_data/NEW_environmental_data/env_layer_names.csv"
    )
  ) %>% tibble() %>% filter(variable_abbreviation %in% names(env))
env_layer_names

# select variable layers from the env rasterstack that are relevant to the study species (edit as needed to exclude any variable names)

# bio_env_layer_names <- env_layer_names

reduced_bio_env_layer_names <- c(
  "AHM",
  "bdod_0_5cm_mean_4km",
  "cfvo_0_5cm_mean_4km",
  "consensus_full_class_1_evergreen_deciduous_needleleaf_trees",
  "consensus_full_class_11_barren",
  "consensus_full_class_2_evergreen_broadleaf_trees",
  "consensus_full_class_3_deciduous_broadleaf_trees",
  "consensus_full_class_4_mixed_other_trees",
  "consensus_full_class_5_shrubs",
  "consensus_full_class_6_herbaceous_vegetation",
  "consensus_full_class_7_cultivated_and_managed_vegetation",
  "consensus_full_class_8_regularly_flooded_vegetation",
  "NFFD_sm",
  "nitrogen_0_5cm_mean_4km",
  "PAS_at",
  "PAS_sm",
  "phh2o_0_5cm_mean_4km",
  "SHM",
  "soc_0_5cm_mean_4km",
  "water_proximity_4km",
  "wv0010_0_5cm_mean_4km",
  "wv0033_0_5cm_mean_4km",
  "wv1500_0_5cm_mean_4km"
)

bio_env_layer_names <-
  env_layer_names %>% filter(variable_abbreviation %in% reduced_bio_env_layer_names)
bio_env_layer_names

# create a subset of the env rasterstack by the selected layer names
bio_env <-
  raster::subset(
    env,
    bio_env_layer_names$variable_abbreviation,
    filename = str_c(raster_tempdir, "bio_env.grd"),
    overwrite = TRUE
  )
names(bio_env)
bio_env

# load elevation raster as a separate layer; this will be used as a reference layer to help visualize the pre-clipped environmental data
elevation <-
  raster::raster(
    str_c(
      main_path,
      "environmental_data/NEW_environmental_data/EarthEnv_Topography/final/elevation_4KMmd_GMTEDmd.tif"
    )
  )
elevation


### Define study area ----

# prepare reference vector layers
boundaries_countries <- read_sf(
  str_c(
    main_path,
    "environmental_data/study_areas/reference_boundaries/Global_Countries.gpkg"
  )
) %>%
  st_transform(st_crs(elevation)) %>%
  st_make_valid() %>%
  dplyr::select(GEOUNIT)
boundaries_states_provinces <- read_sf(
  str_c(
    main_path,
    "environmental_data/study_areas/reference_boundaries/Global_States_Provinces.gpkg"
  )
) %>%
  st_transform(st_crs(elevation)) %>%
  st_make_valid() %>%
  filter(geonunit == "Canada" |
           geonunit == "Mexico" |
           geonunit == "United States of America") %>%
  dplyr::select(name)

# import vector mask layer representing the study area (derived from ecoregions intersecting and surrounding presence locations)
study_area_mask <-
  read_sf(
    "data/study_areas/Ecoregions_Study_Area_culicoides_variipennis_Extended_Mask.gpkg"
  ) %>%
  st_transform(st_crs(elevation)) %>%
  st_make_valid()

# make sf object versions for thinned and unthinned presence coords
presence_sf <- st_as_sf(
  x = presence,
  coords = c("longitude", "latitude"),
  crs = st_crs(elevation)
)
presence_unthinned_sf <- st_as_sf(
  x = presence_unthinned,
  coords = c("longitude", "latitude"),
  crs = st_crs(elevation)
)

# create plot layer from country and state/province boundaries
boundaries_plot_layer <- {
  tm_shape(
    boundaries_countries,
    projection = st_crs(elevation),
    bbox = st_bbox(elevation)
  ) +
    tm_polygons(alpha = 0, border.alpha = 0) +  # transparent copy of polygon (baseline)
    tm_graticules(alpha = 0.2,
                  ticks = FALSE) +
    tm_layout(inner.margins = 0) +  # add graticules (long-lat lines) first
    tm_shape(
      boundaries_countries,
      projection = st_crs(elevation),
      bbox = st_bbox(elevation)
    ) +
    tm_polygons(col = "white", border.alpha = 0) +  # opaque white copy of polygon (to hide graticules under polygon)
    tm_shape(
      boundaries_countries,
      projection = st_crs(elevation),
      bbox = st_bbox(elevation)
    ) +
    tm_polygons(
      alpha = 0.15,
      border.alpha = 1,
      border.col = "black",
      lwd = 1
    ) +  # visible copy of country boundaries
    tm_shape(
      boundaries_states_provinces,
      projection = st_crs(elevation),
      bbox = st_bbox(elevation)
    ) +
    tm_polygons(
      alpha = 0,
      border.alpha = 1,
      border.col = "black",
      lwd = 0.8
    )  # visible copy of state/province boundaries
}

# create plot layer from base fill of elevation raster
elevation_base_plot_layer <- {
  tm_shape(elevation,
           raster.downsample = FALSE) +
    tm_raster(
      palette = "lightgrey",
      legend.show = FALSE,
      alpha = 0.55
    )  # opaque grey base fill of elevation raster
}

blank_map_full_white <- boundaries_plot_layer
blank_map_full_grey <-
  {
    boundaries_plot_layer + elevation_base_plot_layer
  }

# save the image output
tmap_save(blank_map_full_white,
          filename = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/blank_map_full_white.png",
          dpi = 400)

# save the image output
tmap_save(blank_map_full_grey,
          filename = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/blank_map_full_grey.png",
          dpi = 400)

# save the image output
tmap_save(blank_map_full_white,
          filename = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/blank_map_full_white.pdf")

# save the image output
tmap_save(blank_map_full_grey,
          filename = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/blank_map_full_grey.pdf")

# map out the presence locations prior to establishing the study area
(presence_unclipped_plot <- {
  boundaries_plot_layer +
    elevation_base_plot_layer +
    tm_shape(presence_sf) +
    tm_dots(
      size = 0.08,
      shape = 21,
      col = "red",
      alpha = 0.4,
      border.lwd = 0.5,
      border.alpha = 0.9
    )
})

# save the image output
tmap_save(presence_unclipped_plot,
          filename = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/presence_unclipped_plot.png",
          dpi = 400)

# create plot layer from study area mask
study_area_mask_plot_layer <- {
  tm_shape(
    study_area_mask,
    projection = st_crs(study_area_mask),
    bbox = study_area_mask,
    is.master = TRUE
  ) +
    tm_polygons(
      col = "cornsilk",
      alpha = 0.3,
      border.alpha = 0.5,
      lwd = 1.8
    )
}

# plot the presence locations within the study area
(study_area_plot <- {
  boundaries_plot_layer +
    elevation_base_plot_layer +
    study_area_mask_plot_layer +
    tm_shape(presence_sf) +
    tm_dots(
      size = 0.15,
      shape = 21,
      col = "red",
      alpha = 0.4,
      border.lwd = 0.5,
      border.alpha = 0.9
    )
})

# save the image output
tmap_save(study_area_plot,
          filename = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/study_area_plot.png",
          dpi = 400)

# clip the environmental rasterstack and elevation layer by the study area mask layer
bio_env_c <-
  raster::mask(bio_env, study_area_mask) %>% raster::trim(filename = str_c(raster_tempdir, "bio_env_c.grd"),
                                                          overwrite = TRUE)
elevation_c <-
  raster::mask(elevation, study_area_mask) %>% raster::trim(filename = str_c(raster_tempdir, "elevation_c.grd"),
                                                            overwrite = TRUE)

# create plot layer from base fill of clipped elevation raster
elevation_c_base_plot_layer <- {
  tm_shape(
    elevation_c,
    raster.downsample = FALSE,
    projection = st_crs(study_area_mask),
    bbox = study_area_mask,
    is.master = TRUE
  ) +
    tm_raster(
      palette = "lightgrey",
      legend.show = FALSE,
      alpha = 0.55
    )  # opaque grey base fill of elevation raster
}

# plot the presence locations with the background reference raster clipped to the study area
(presence_clipped_plot <- {
  boundaries_plot_layer +
    elevation_c_base_plot_layer +
    tm_shape(presence_sf) +
    tm_dots(
      size = 0.15,
      shape = 21,
      col = "red",
      alpha = 0.4,
      border.lwd = 0.5,
      border.alpha = 0.9
    )
})

# save the image output
tmap_save(presence_clipped_plot,
          filename = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/presence_clipped_plot.png",
          dpi = 400)

# plot a map showing the clipped elevation layer together with the presence locations, including points from both before and after thinning for comparison
(presence_thinning_comparison_plot <- {
  boundaries_plot_layer +
    elevation_c_base_plot_layer +
    tm_shape(presence_unthinned_sf) +
    tm_dots(
      size = 0.15,
      shape = 21,
      col = "royalblue4",
      border.lwd = 0.3,
      border.col = "navyblue"
    ) +
    tm_shape(presence_sf) +
    tm_dots(
      size = 0.15,
      shape = 21,
      col = "indianred3",
      border.lwd = 0.3,
      border.col = "darkred"
    ) +
    tm_add_legend(
      type = "symbol",
      labels = c("After spatial filtering", "Before spatial filtering"),
      shape = 21,
      col = c("indianred3", "royalblue4"),
      border.lwd = 0.3,
      border.col = c("darkred", "navyblue")
    ) +
    tm_layout(legend.position = c("left", "bottom"))
})

# save the image output
tmap_save(presence_thinning_comparison_plot,
          filename = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/presence_thinning_comparison_plot.png",
          dpi = 400)


### Generate background locations ----

# determine number of random background points to be generated per decade based on the corresponding decade's proportion of collected presence locations
bg_sample_size_by_timestep <-
  floor(10000 * data.frame(prob_table_timestep)$Freq)

# through independent random sampling, generate background data corresponding to each record collection decade, where the number of background points per set is roughly proportional to the presence record frequency for each decade
set.seed(25)
bg_list <-
  lapply(bg_sample_size_by_timestep,
         randomPoints,
         mask = bio_env_c,
         extf = 1)

# rename time-specific background data subsets in list for clarity
names(bg_list) <- paste0("bg_", names(prob_table_timestep))
summary(bg_list)

# save all objects in R environment
save.image(
  "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil.RData"
)


### Obtain time-specific presence locations and environmental data ----

# define timesteps as vector object
timesteps <- prob_table_timestep %>% names

# build custom function to group time-specific presence localities and environmental data
get.presence.and.env.by.timestep <- function(timesteps) {
  timestep <- timesteps

  ### Prepare presence locations ----

  # filter the presence locations by a given collection timestep
  presence_data_filtered_by_timestep <-
    presence_data_filtered %>% filter(decade == timestep)

  # get just coordinate values from occurrence data set
  presence_by_timestep <-
    presence_data_filtered_by_timestep %>% dplyr::select(longitude, latitude)

  ### Acquire environmental variables ----

  # import raster files (environmental variables/predictors)
  files_by_timestep <-
    c(
      list.files(
        dir(
          str_c(
            main_path,
            "environmental_data/NEW_environmental_data/ClimateNA/final/decadal_1911_2020"
          ),
          pattern = str_c("Decade_", str_sub(timestep, end = 3)),
          full.names = TRUE,
          recursive = TRUE,
          include.dirs = TRUE
        ),
        pattern = "tif$",
        full.names = TRUE,
        recursive = TRUE
      ),
      list.files(
        c(
          str_c(
            main_path,
            "environmental_data/NEW_environmental_data/EarthEnv_Consensus_Land_Cover/final"
          ),
          str_c(
            main_path,
            "environmental_data/NEW_environmental_data/HydroSHEDS_Hydrology/final"
          ),
          str_c(
            main_path,
            "environmental_data/NEW_environmental_data/SoilGrids/final_NAfilled"
          )
        ),
        pattern = 'tif$',
        full.names = TRUE,
        recursive = TRUE
      )
    )
  files_by_timestep <-
    setdiff(files_by_timestep,
            str_subset(files_by_timestep, "monthly"))

  ### Define study area ----

  # convert the imported files into a rasterstack object, create a subset of it by the predefined selected layer names, and then clip it by the vector mask layer
  bio_env_c_by_timestep <-
    files_by_timestep %>%
    raster::stack() %>%
    raster::subset(bio_env_layer_names$variable_abbreviation) %>%
    raster::mask(study_area_mask) %>%
    raster::trim(
      filename = str_c(raster_tempdir, "bio_env_c_", timestep, ".grd"),
      overwrite = TRUE
    )

  # make sf object version for presence coords by timestep
  presence_by_timestep_sf <- st_as_sf(
    x = presence_by_timestep,
    coords = c("longitude", "latitude"),
    crs = st_crs(study_area_mask)
  )

  # create presence by timestep plot layer
  presence_by_timestep_plot_layer <- {
    tm_shape(
      presence_by_timestep_sf,
      projection = st_crs(study_area_mask),
      bbox = study_area_mask,
      is.master = TRUE
    ) +
      tm_dots(
        size = 0.15,
        shape = 21,
        col = "red",
        alpha = 0.4,
        border.lwd = 0.5,
        border.alpha = 0.9
      ) +  # presences points
      tm_layout(
        panel.labels = str_c(
          "Presence localities (",
          timestep %>% as.numeric %>% +1 %>% as.character,
          "-",
          timestep %>% as.numeric %>% +10 %>% as.character,
          ")"
        )
      )  # plot title (panel label format)
  }

  # plot presence localities by timestep
  presence_by_timestep_plot <-
    boundaries_plot_layer +
    elevation_c_base_plot_layer +
    presence_by_timestep_plot_layer

  # save the image output
  tmap_save(
    presence_by_timestep_plot,
    filename = str_c(
      "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/presence_by_timestep_plot_",
      timestep,
      ".png"
    ),
    dpi = 400
  )

  # return list of variables to pass through to output object
  list(
    presence_data_filtered_by_timestep = presence_data_filtered_by_timestep,
    presence_by_timestep = presence_by_timestep,
    files_by_timestep = files_by_timestep,
    bio_env_c_by_timestep = bio_env_c_by_timestep
  )
}

# extract time-specific environmental data at presence localities using custom function

{
  get_presence_and_env_by_timestep_output_start_time <- proc.time()

  {
    plan(multisession, workers = num_parallel_workers)
  }  # initiate parallel processing

  get_presence_and_env_by_timestep_output <-
    future_sapply(
      timesteps,
      get.presence.and.env.by.timestep,
      simplify = FALSE,
      USE.NAMES = TRUE
    ) %>% list_transpose()

  {
    plan(sequential)
    }  # end parallel processing

  get_presence_and_env_by_timestep_output_stop_time <-
    proc.time()
  get_presence_and_env_by_timestep_output_run_time <-
    hms::hms(
      get_presence_and_env_by_timestep_output_stop_time[["elapsed"]] - get_presence_and_env_by_timestep_output_start_time[["elapsed"]]
    )
  print(str_c(
    "Run time (HMS): ",
    get_presence_and_env_by_timestep_output_run_time
  ))
}

summary(get_presence_and_env_by_timestep_output)

# save object to disk
save(get_presence_and_env_by_timestep_output, file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/get_presence_and_env_by_timestep_output.RData")

# extract time-specific presence coords from output and store them in new list, then rename elements for clarity
presence_list <-
  get_presence_and_env_by_timestep_output$presence_by_timestep
names(presence_list) <- paste0("presence_", names(presence_list))
summary(presence_list)

# extract time-specific environmental data from output and store them in new list, then rename elements for clarity
bio_env_c_list <-
  get_presence_and_env_by_timestep_output$bio_env_c_by_timestep
names(bio_env_c_list) <- paste0("bio_env_c_", names(bio_env_c_list))
summary(bio_env_c_list)

# save all objects in R environment
save.image(
  "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil.RData"
)

### Create SWD object ----

# use the prepareSWD function to create an SWD object that stores the species name, coordinates of presence and background locations, and value of environmental variables at those locations (excluding those that have NA values for at least one environmental variable)

{
  devtools::dev_mode()
  {
    p_unload(SDMtune)
    p_unlock()
    if (p_isinstalled(SDMtune) == FALSE |
        (p_isinstalled(SDMtune) == TRUE &&
         p_version(SDMtune) != "1.1.6")) {
      devtools::install_version("SDMtune", "1.1.6")
    }
    p_unlock()
    library(SDMtune, lib.loc = .libPaths()[1])
  }
}  # unload current SDMtune package and instead use older version (1.1.6) to prevent R crashing while running next command (due to GDAL error)

data_list_before_addSamplesToBg <- mapply(
  prepareSWD,
  species = "Culicoides variipennis",
  p = presence_list,
  a = bg_list,
  env = bio_env_c_list,
  SIMPLIFY = FALSE,
  USE.NAMES = FALSE
)

{
  p_unload(SDMtune)
  devtools::dev_mode()
  p_load(SDMtune)
}  # revert back to SDMtune package version 1.2.0

# add presence data to background
data_list <- lapply(data_list_before_addSamplesToBg, addSamplesToBg)

# extract individual time-specific SWD object subsets from the list
names(data_list) <- paste0("data_", names(prob_table_timestep))
summary(data_list)

# merge the SWD objects together; this should result in time-specific presence, background, and environmental data all pooled together into one large SWD object
data <- mergeSWD(data_list[[1]],
                 data_list[[2]])
for (i in seq_along(data_list)) {
  if (i >= 3)
  {
    data <- mergeSWD(data,
                     data_list[[i]])
  }
}

# check out the created SWD object
data
summary(data)

# check the environmental values
head(data@data)

# check the coordinates
head(data@coords)

# check the species label
data@species

# save the object in a single file with the column pa indicating if the location is a presence (1) or an absence/background (0) site
swd2csv(data, file_name = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/data.csv")

# save all objects in R environment
save.image(
  "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil.RData"
)


# Establish nested cross-validation structure for model evaluation ----

## Build outer loop cross-validation folds to split data into train and test sets ----

# define number of cross-validation folds in outer loop
num_outer_k_folds <- 5

# create the outer loop cross-validation folds
outerloop_cv_folds <- randomFolds(data,
                                  k = num_outer_k_folds,
                                  only_presence = TRUE,
                                  seed = 100)
str(outerloop_cv_folds)

# save train and test outer folds into separate matrices
zeallot::`%<-%`(c(outerloop_cv_train_folds, outerloop_cv_test_folds),
                outerloop_cv_folds)

# create function to group full SWD object by outer folds; the goal is to create 5 subsetted SWD objects total, each excluding a different test fold
subset_swd_by_outerfolds <- function(swd, fold) {
  data <- swd@data[fold, , drop = FALSE]
  coords <- swd@coords[fold, , drop = FALSE]
  rownames(data) <- NULL
  rownames(coords) <- NULL
  pa <- swd@pa[fold]

  SWD(
    species = swd@species,
    data = data,
    coords = coords,
    pa = pa
  )
}

# run subset_swd_by_outerfolds function on training data for each set of outer loop cv folds to obtain the 5 subsetted SWD objects in the outer loop
outerloop_train_subset_list <-
  sapply(1:num_outer_k_folds, function(x)
    NULL)

for (j in 1:num_outer_k_folds) {
  output <-
    subset_swd_by_outerfolds(swd = data, fold = outerloop_cv_train_folds[, j])
  outerloop_train_subset_list[[j]] <- output
}
outerloop_train_subset_list
summary(outerloop_train_subset_list)

# run subset_swd_by_outerfolds function on test data for each set of outer loop cv folds to obtain the 5 subsetted SWD objects in the outer loop
outerloop_test_subset_list <-
  sapply(1:num_outer_k_folds, function(x)
    NULL)

for (j in 1:num_outer_k_folds) {
  output <-
    subset_swd_by_outerfolds(swd = data, fold = outerloop_cv_test_folds[, j])
  outerloop_test_subset_list[[j]] <- output
}

outerloop_test_subset_list
summary(outerloop_test_subset_list)

# save the train and test data sets as csv files
for (i in 1:num_outer_k_folds) {
  swd2csv(
    outerloop_train_subset_list[[i]],
    file_name = str_c(
      "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/outerloop_train_subset_",
      i,
      ".csv"
    )
  )
  swd2csv(
    outerloop_test_subset_list[[i]],
    file_name = str_c(
      "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/outerloop_test_subset_",
      i,
      ".csv"
    )
  )
}

# save all objects in R environment
save.image(
  "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil.RData"
)

# visualize the outer loop train and test subsets (where partitions 1 and 2 represent train and test, respectively)
train_test_subset_plot_list <-
  sapply(1:num_outer_k_folds, function(x)
    NULL)

{
  cl <- makeCluster(num_parallel_workers)
  registerDoParallel(cl)
  }  # initiate parallel processing

train_test_subset_plot_list <-
  foreach(i = 1:num_outer_k_folds,
          .packages = c("ENMeval", "tidyverse")) %dopar% {
            train_test_subset_plot <-
              evalplot.grps(
                pts = rbind(
                  outerloop_train_subset_list[[i]]@coords[outerloop_train_subset_list[[i]]@pa == 1, ],
                  outerloop_test_subset_list[[i]]@coords[outerloop_test_subset_list[[i]]@pa == 1, ]
                ),
                pts.grp = c(rep(
                  1, nrow(outerloop_train_subset_list[[i]]@coords[outerloop_train_subset_list[[i]]@pa == 1, ])
                ), rep(
                  2, nrow(outerloop_test_subset_list[[i]]@coords[outerloop_test_subset_list[[i]]@pa == 1, ])
                )),
                envs = elevation_c,
                pts.size = 0.5
              )

            # save the image output
            ggsave(
              file = str_c("train_test_subset_plot_", i, ".png"),
              train_test_subset_plot,
              device = "png",
              path = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil",
              width = 7,
              height = 7,
              units = "in",
              dpi = 400
            )
            ggsave(
              file = str_c("train_test_subset_plot_", i, ".pdf"),
              train_test_subset_plot,
              device = "pdf",
              path = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil",
              width = 7,
              height = 7,
              units = "in"
            )

            train_test_subset_plot
          }

{
  stopCluster(cl)
  registerDoSEQ()
}  # end parallel processing

# save all objects in R environment
save.image(
  "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil.RData"
)


## Build inner loop cross-validation folds to split each outer loop train data set into train and validation sets ----

# for each outer loop cross-validation group, run inner loop cross-validation (k = 4) using training data sets derived from outer loop cross-validation folds (k = 5), then get outer loop cross-validation output
innerloop_cv_folds_list <- sapply(1:num_outer_k_folds, function(x)
  NULL)
innerloop_cv_model_untuned_list <-
  sapply(1:num_outer_k_folds, function(x)
    NULL)

{
  cl <- makeCluster(num_parallel_workers)
  registerDoParallel(cl)
  }  # initiate parallel processing

innerloop_cv_folds_list <-
  foreach(i = 1:num_outer_k_folds,
          .packages = c("ENMeval")) %dopar% {
            # create cross-validation folds (k = 4) for training data set
            set.seed(25)
            innerloop_cv_folds <-
              get.checkerboard2(
                occs = outerloop_train_subset_list[[i]]@coords[outerloop_train_subset_list[[i]]@pa == 1,],
                envs = bio_env_c,
                bg = outerloop_train_subset_list[[i]]@coords[outerloop_train_subset_list[[i]]@pa == 0,],
                aggregation.factor = 4
              )
            innerloop_cv_folds
          }
innerloop_cv_folds_list
summary(innerloop_cv_folds_list)

if (Sys.info()[["sysname"]] == "Windows") {
  set.seed(25)
  innerloop_cv_model_untuned_list <- mapply(train,
                                            data = outerloop_train_subset_list,
                                            method = "Maxent",
                                            folds = innerloop_cv_folds_list)
  innerloop_cv_model_untuned_list
  summary(innerloop_cv_model_untuned_list)
} else {
  innerloop_cv_model_untuned_list <-
    foreach(i = 1:num_outer_k_folds,
            .packages = c("SDMtune")) %dopar% {
              # train model via cross-validation
              set.seed(25)
              innerloop_cv_model_untuned <-
                train(outerloop_train_subset_list[[i]],
                      method = "Maxent",
                      folds = innerloop_cv_folds_list[[i]])
              innerloop_cv_model_untuned
            }
  innerloop_cv_model_untuned_list
  summary(innerloop_cv_model_untuned_list)
}

{
  stopCluster(cl)
  registerDoSEQ()
}  # end parallel processing

# save all objects in R environment
save.image(
  "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil.RData"
)

# plot the inner loop checkerboard partitions (k = 4) of the presence points
innerloop_cv_folds_presence_plot_list <-
  sapply(1:num_outer_k_folds, function(x)
    NULL)

{
  cl <- makeCluster(num_parallel_workers)
  registerDoParallel(cl)
  }  # initiate parallel processing

innerloop_cv_folds_presence_plot_list <-
  foreach(i = 1:num_outer_k_folds,
          .packages = c("ENMeval", "tidyverse")) %dopar% {
            innerloop_cv_folds_presence_plot <-
              evalplot.grps(
                pts = outerloop_train_subset_list[[i]]@coords[outerloop_train_subset_list[[i]]@pa == 1,],
                pts.grp = innerloop_cv_folds_list[[i]]$occs.grp,
                envs = elevation_c,
                pts.size = 0.5
              )

            # save the image output
            ggsave(
              file = str_c("innerloop_cv_folds_presence_plot_", i, ".png"),
              innerloop_cv_folds_presence_plot,
              device = "png",
              path = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil",
              width = 7,
              height = 7,
              units = "in",
              dpi = 400
            )
            ggsave(
              file = str_c("innerloop_cv_folds_presence_plot_", i, ".pdf"),
              innerloop_cv_folds_presence_plot,
              device = "pdf",
              path = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil",
              width = 7,
              height = 7,
              units = "in"
            )

            innerloop_cv_folds_presence_plot
          }

{
  stopCluster(cl)
  registerDoSEQ()
}  # end parallel processing

# save all objects in R environment
save.image(
  "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil.RData"
)

# plot the inner loop checkerboard partitions (k = 4) of the background points
innerloop_cv_folds_background_plot_list <-
  sapply(1:num_outer_k_folds, function(x)
    NULL)

{
  cl <- makeCluster(num_parallel_workers)
  registerDoParallel(cl)
  }  # initiate parallel processing

innerloop_cv_folds_background_plot_list <-
  foreach(i = 1:num_outer_k_folds,
          .packages = c("ENMeval", "tidyverse")) %dopar% {
            innerloop_cv_folds_background_plot <-
              evalplot.grps(
                pts = outerloop_train_subset_list[[i]]@coords[outerloop_train_subset_list[[i]]@pa == 0, ],
                pts.grp = innerloop_cv_folds_list[[i]]$bg.grp,
                envs = elevation_c,
                pts.size = 0.5
              )

            # save the image output
            ggsave(
              file = str_c("innerloop_cv_folds_background_plot_", i, ".png"),
              innerloop_cv_folds_background_plot,
              device = "png",
              path = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil",
              width = 7,
              height = 7,
              units = "in",
              dpi = 400
            )
            ggsave(
              file = str_c("innerloop_cv_folds_background_plot_", i, ".pdf"),
              innerloop_cv_folds_background_plot,
              device = "pdf",
              path = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil",
              width = 7,
              height = 7,
              units = "in"
            )

            innerloop_cv_folds_background_plot
          }

{
  stopCluster(cl)
  registerDoSEQ()
}  # end parallel processing

# save all objects in R environment
save.image(
  "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil.RData"
)


## Tune model hyperparameters ----

# see which hyperparameters can be tuned based on the given model's algorithm
getTunableArgs(innerloop_cv_model_untuned_list[[1]])

# create a list of hyperparameter combinations to be tested, including L1 regularization multiplier values, feature classes, and number of iterations
hyperparam_combos <-
  list(
    reg = seq(0.1, 10, 0.1),
    fc = c("l", "lq", "lh", "lp", "lqp", "lqph"),
    iter = seq(500, 2000, 100)
  )
summary(hyperparam_combos)


### Optimize search (genetic algorithm) ----

# perform an optimize search to find the best possible combination of hyperparameters for the model chosen via a data-driven genetic algorithm from an initial random sample
innerloop_cv_optimizesearch_output_list <-
  sapply(1:num_outer_k_folds, function(x)
    NULL)

{
  innerloop_cv_optimizesearch_output_list_start_time <- proc.time()

  if (Sys.info()[["sysname"]] == "Windows") {
    innerloop_cv_optimizesearch_output_list <- lapply(
      innerloop_cv_model_untuned_list,
      optimizeModel,
      hypers = hyperparam_combos,
      metric = "auc",
      pop = 30,
      gen = 5,
      interactive = FALSE,
      seed = 25
    )
  } else {
    {
      cl <- makeCluster(num_parallel_workers)
      registerDoParallel(cl)
    }  # initiate parallel processing

    innerloop_cv_optimizesearch_output_list <-
      parLapply(
        cl,
        innerloop_cv_model_untuned_list,
        optimizeModel,
        hypers = hyperparam_combos,
        metric = "auc",
        pop = 30,
        gen = 5,
        interactive = FALSE,
        seed = 25
      )

    {
      stopCluster(cl)
      registerDoSEQ()
      }  # end parallel processing
  }

  innerloop_cv_optimizesearch_output_list_stop_time <-
    proc.time()
  innerloop_cv_optimizesearch_output_list_run_time <-
    hms::hms(
      innerloop_cv_optimizesearch_output_list_stop_time[["elapsed"]] - innerloop_cv_optimizesearch_output_list_start_time[["elapsed"]]
    )
  print(str_c(
    "Run time (HMS): ",
    innerloop_cv_optimizesearch_output_list_run_time
  ))
  }

innerloop_cv_optimizesearch_output_list
summary(innerloop_cv_optimizesearch_output_list)

# save object to disk
save(innerloop_cv_optimizesearch_output_list, file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/innerloop_cv_optimizesearch_output_list.RData")

# save all objects in R environment
save.image(
  "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil.RData"
)

innerloop_cv_optimizesearch_output_results_list <-
  sapply(1:num_outer_k_folds, function(x)
    NULL)
innerloop_cv_optimizesearch_output_plot_list <-
  sapply(1:num_outer_k_folds, function(x)
    NULL)
hyperparam_optimizesearch_list <-
  sapply(1:num_outer_k_folds, function(x)
    NULL)
innerloop_cv_model_tuned_list <-
  sapply(1:num_outer_k_folds, function(x)
    NULL)
innerloop_cv_model_tuned_data_list <-
  sapply(1:num_outer_k_folds, function(x)
    NULL)

for (i in 1:num_outer_k_folds) {
  # check the output (ordered by test AUC in descending order by default), including only unique results
  (
    innerloop_cv_optimizesearch_output_results_list[[i]] <-
      innerloop_cv_optimizesearch_output_list[[i]]@results %>% unique()
  )

  # plot the hyperparameter tuning search chart
  innerloop_cv_optimizesearch_output_plot_list[[i]] <-
    plot(innerloop_cv_optimizesearch_output_list[[i]],
         title = "Hyperparameter Tuning:\nOptimize Search")

  # save the image output
  ggsave(
    file = str_c("innerloop_cv_optimizesearch_output_plot_",
                 i,
                 ".pdf"),
    innerloop_cv_optimizesearch_output_plot_list[[i]],
    device = "pdf",
    path = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil",
    width = 7,
    height = 6,
    units = "in"
  )
  ggsave(
    file = str_c("innerloop_cv_optimizesearch_output_plot_",
                 i,
                 ".png"),
    innerloop_cv_optimizesearch_output_plot_list[[i]],
    device = "png",
    path = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil",
    width = 7,
    height = 6,
    units = "in",
    dpi = 400
  )

  # save the best hyperparameter combination to an object
  (
    hyperparam_optimizesearch_list[[i]] <-
      innerloop_cv_optimizesearch_output_results_list[[i]][1,]
  )

  # save the best hyperparameter combination as a csv file
  write.csv(
    hyperparam_optimizesearch_list[[i]],
    file = str_c(
      "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/hyperparam_optimizesearch_",
      i,
      ".csv"
    ),
    row.names = FALSE
  )

  # pull the best cv model from the tuning results and assign it to new model object in list
  (
    innerloop_cv_model_tuned_list[[i]] <-
      innerloop_cv_optimizesearch_output_list[[i]]@models[[i]]
  )

  # save the output from the new tuned cv model as a txt file
  sink(
    str_c(
      "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/innerloop_cv_model_tuned_",
      i,
      ".txt"
    )
  )
  print(innerloop_cv_model_tuned_list[[i]])
  sink()

  innerloop_cv_model_tuned_data_list[[i]] <-
    innerloop_cv_model_tuned_list[[i]]@data

  # save data from the new tuned cv model as a csv file
  swd2csv(
    innerloop_cv_model_tuned_data_list[[i]],
    file_name = str_c(
      "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/innerloop_cv_model_tuned_data_",
      i,
      ".csv"
    )
  )
}

# save all objects in R environment
save.image(
  "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil.RData"
)


## Evaluate inner loop cross-validation procedure using withheld outer loop test data ----

# build each tuned model using entire corresponding outer loop train data set (without cross-validation)
outerloop_model_tuned_list <-
  sapply(1:num_outer_k_folds, function(x)
    NULL)

outerloop_model_tuned_list <- foreach(i = 1:num_outer_k_folds,
                                      .packages = "SDMtune") %do% {
                                        set.seed(25)
                                        train(
                                          data = innerloop_cv_model_tuned_data_list[[i]],
                                          method = "Maxent",
                                          fc = hyperparam_optimizesearch_list[[i]][, "fc"],
                                          reg = hyperparam_optimizesearch_list[[i]][, "reg"],
                                          iter = hyperparam_optimizesearch_list[[i]][, "iter"]
                                        )
                                      }

outerloop_model_tuned_list

# for each tuned outer loop model, create ROC (Receiver Operator Characteristic) plot showing train and test AUC, using corresponding outer loop test data set to evaluate the model and calculate the test AUC
outerloop_cv_output_stats_ROC_plot_list <-
  sapply(1:num_outer_k_folds, function(x)
    NULL)

for (i in 1:num_outer_k_folds) {
  outerloop_cv_output_stats_ROC_plot_list[[i]] <-
    plotROC(outerloop_model_tuned_list[[i]], test = outerloop_test_subset_list[[i]])

  # save the image output
  ggsave(
    file = str_c("outerloop_cv_output_stats_ROC_plot_",
                 i,
                 ".pdf"),
    outerloop_cv_output_stats_ROC_plot_list[[i]],
    device = "pdf",
    path = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil",
    width = 7,
    height = 6,
    units = "in"
  )
  ggsave(
    file = str_c("outerloop_cv_output_stats_ROC_plot_",
                 i,
                 ".png"),
    outerloop_cv_output_stats_ROC_plot_list[[i]],
    device = "png",
    path = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil",
    width = 7,
    height = 6,
    units = "in",
    dpi = 400
  )
}

# create empty data frame for evaluation output from outer loop cross-validation with set number of folds
outerloop_cv_output_stats <-
  data.frame(matrix(NA, nrow = num_outer_k_folds, ncol = 14))

colnames(outerloop_cv_output_stats) <-
  c(
    "tuned_model_id",
    "fc",
    "reg",
    "iter",
    "AUC",
    "AUC_diff",
    "AUC_diff_value_sign",
    "TSS",
    "TSS_diff",
    "TSS_diff_value_sign",
    "OR_min_train_presence",
    "OR_max_train_se_plus_sp",
    "OR_diff_max_train_se_plus_sp",
    "OR_diff_max_train_se_plus_sp_value_sign"
  )

outerloop_cv_output_stats
str(outerloop_cv_output_stats)

# calculate outer loop evaluation metrics using withheld test data to determine performance of cross-validation model tuning (k = 5)
for (i in 1:num_outer_k_folds) {
  # get number/ID of tuned model
  outerloop_cv_output_stats$tuned_model_id[i] <-
    as.character(i)

  # get optimal hyperparameter setting for feature classes
  outerloop_cv_output_stats$fc[i] <-
    hyperparam_optimizesearch_list[[i]][, "fc"]

  # get optimal hyperparameter setting for regularization multiplier
  outerloop_cv_output_stats$reg[i] <-
    as.character(hyperparam_optimizesearch_list[[i]][, "reg"])

  # get optimal hyperparameter setting for maximum number of iterations
  outerloop_cv_output_stats$iter[i] <-
    as.character(hyperparam_optimizesearch_list[[i]][, "iter"])

  # calculate test AUC
  outerloop_cv_output_stats$AUC[i] <-
    auc(outerloop_model_tuned_list[[i]], test = outerloop_test_subset_list[[i]])

  # calculate absolute difference between train and test AUC
  outerloop_cv_output_stats$AUC_diff[i] <-
    abs(
      auc(outerloop_model_tuned_list[[i]]) - auc(outerloop_model_tuned_list[[i]], test = outerloop_test_subset_list[[i]])
    )

  # determine value sign of difference between train and test AUC; options include positive (model overfit to test data), negative (model underfit to test data), or zero (model perfectly fit to test data)
  outerloop_cv_output_stats$AUC_diff_value_sign[i] <-
    switch(
      as.character(sign(
        auc(outerloop_model_tuned_list[[i]]) - auc(outerloop_model_tuned_list[[i]], test = outerloop_test_subset_list[[i]])
      )),
      "1" = "Positive (+)",
      "-1" = "Negative (-)",
      "0" = "Zero (0)",
      NULL
    )

  # calculate test TSS
  outerloop_cv_output_stats$TSS[i] <-
    tss(outerloop_model_tuned_list[[i]], test = outerloop_test_subset_list[[i]])

  # calculate absolute difference between train and test TSS
  outerloop_cv_output_stats$TSS_diff[i] <-
    abs(
      tss(outerloop_model_tuned_list[[i]]) - tss(outerloop_model_tuned_list[[i]], test = outerloop_test_subset_list[[i]])
    )

  # determine value sign of difference between train and test TSS; options include positive (model overfit to test data), negative (model underfit to test data), or zero (model perfectly fit to test data)
  outerloop_cv_output_stats$TSS_diff_value_sign[i] <-
    switch(
      as.character(sign(
        tss(outerloop_model_tuned_list[[i]]) - tss(outerloop_model_tuned_list[[i]], test = outerloop_test_subset_list[[i]])
      )),
      "1" = "Positive (+)",
      "-1" = "Negative (-)",
      "0" = "Zero (0)",
      NULL
    )

  # calculate test omission rate at minimum training presence threshold
  outerloop_cv_output_stats$OR_min_train_presence[i] <-
    thresholds(outerloop_model_tuned_list[[i]],
               type = "cloglog",
               test = outerloop_test_subset_list[[i]])[1, 5]

  # calculate test omission rate at maximum training sensitivity plus specificity threshold
  outerloop_cv_output_stats$OR_max_train_se_plus_sp[i] <-
    thresholds(outerloop_model_tuned_list[[i]],
               type = "cloglog",
               test = outerloop_test_subset_list[[i]])[3, 5]

  # calculate absolute difference between train and test omission rates at maximum training sensitivity plus specificity threshold
  outerloop_cv_output_stats$OR_diff_max_train_se_plus_sp[i] <-
    abs(
      thresholds(
        outerloop_model_tuned_list[[i]],
        type = "cloglog",
        test = outerloop_test_subset_list[[i]]
      )[3, 4]
      -
        thresholds(
          outerloop_model_tuned_list[[i]],
          type = "cloglog",
          test = outerloop_test_subset_list[[i]]
        )[3, 5]
    )

  # determine value sign of difference between train and test omission rates at maximum training sensitivity plus specificity threshold; options include positive (model overfit to test data), negative (model underfit to test data), or zero (model perfectly fit to test data)
  outerloop_cv_output_stats$OR_diff_max_train_se_plus_sp_value_sign[i] <-
    switch(
      as.character(sign(
        thresholds(
          outerloop_model_tuned_list[[i]],
          type = "cloglog",
          test = outerloop_test_subset_list[[i]]
        )[3, 4]
        -
          thresholds(
            outerloop_model_tuned_list[[i]],
            type = "cloglog",
            test = outerloop_test_subset_list[[i]]
          )[3, 5]
      )),
      "1" = "Positive (+)",
      "-1" = "Negative (-)",
      "0" = "Zero (0)",
      NULL
    )
}

outerloop_cv_output_stats <-
  as_tibble(outerloop_cv_output_stats)

# check the per-model outer loop cv output statistics
outerloop_cv_output_stats
str(outerloop_cv_output_stats)

# save the per-model outer loop cv output statistics as a csv file
write.csv(outerloop_cv_output_stats,
          file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/outerloop_cv_output_stats.csv",
          row.names = FALSE)

# make version of table with values rounded to 3 digits
outerloop_cv_output_stats_rounded <-
  outerloop_cv_output_stats %>%
  mutate(across(where(is.character), ~ as.factor(.x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x >= 0.001, round(.x, 3), .x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x < 0.001 & .x != 0, "<0.001", .x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x == 0, "0", .x)))

# save output as a csv file
write.csv(outerloop_cv_output_stats_rounded,
          file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/outerloop_cv_output_stats_rounded.csv",
          row.names = FALSE)

# calculating summary statistics from outer loop cross-validation output (where sample size (n) = number of outer loop cross-validation folds (k)); these values, particularly the mean and standard deviation, can then be reported alongside the final model as estimated generalization error (estimated performance of final model on new/unseen data)
outerloop_cv_output_stats_summary <-
  outerloop_cv_output_stats %>%
  dplyr::select(where(is.numeric)) %>%  # picking just the numeric-value columns
  pivot_longer(everything(),
               names_to = "performance_metric") %>%  # pivoting data frame (making it long) as prep for calculating descriptive stats
  group_by(performance_metric) %>%
  summarize(
    mean = mean(value) %>% signif(4),
    sd = sd(value) %>% signif(4),
    min = min(value) %>% signif(4),
    max = max(value) %>% signif(4),
    q1 = quantile(value, 0.25) %>% signif(4),
    median = median(value) %>% signif(4),
    q3 = quantile(value, 0.75) %>% signif(4),
    IQR = IQR(value) %>% signif(4)
  ) %>%
  arrange(factor(
    performance_metric,
    levels = c(
      "AUC",
      "AUC_diff",
      "TSS",
      "TSS_diff",
      "OR_min_train_presence",
      "OR_max_train_se_plus_sp",
      "OR_diff_max_train_se_plus_sp"
    )
  ))

# check the outer loop cv output statistical summary (estimated generalization error)
outerloop_cv_output_stats_summary

# save the outer loop cv output statistical summary as a csv file
write.csv(outerloop_cv_output_stats_summary,
          file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/outerloop_cv_output_stats_summary_(estimated_generalization_error).csv",
          row.names = FALSE)

# make version of table with values rounded to 3 digits
outerloop_cv_output_stats_summary_rounded <-
  outerloop_cv_output_stats_summary %>%
  mutate(across(where(is.character), ~ as.factor(.x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x >= 0.001, round(.x, 3), .x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x < 0.001 & .x != 0, "<0.001", .x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x == 0, "0", .x)))
outerloop_cv_output_stats_summary_rounded

# save output as a csv file
write.csv(outerloop_cv_output_stats_summary_rounded,
          file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/outerloop_cv_output_stats_summary_(estimated_generalization_error)_rounded.csv",
          row.names = FALSE)

# make subsetted and cleaned version of rounded table; save output as a csv file
outerloop_cv_output_stats_summary_rounded_main <-
  outerloop_cv_output_stats_summary_rounded %>%
  dplyr::select(performance_metric,
                mean,
                sd,
                min,
                max) %>%
  rename(
    "Performance metric" = "performance_metric",
    "Mean" = "mean",
    "Standard deviation" = "sd",
    "Minimum" = "min",
    "Maximum" = "max"
  )

write.csv(
  outerloop_cv_output_stats_summary_rounded_main,
  file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/outerloop_cv_output_stats_summary_(estimated_generalization_error)_rounded_main.csv",
  row.names = FALSE
)

outerloop_cv_output_stats_summary_rounded_main_clean <-
  outerloop_cv_output_stats_summary_rounded_main

outerloop_cv_output_stats_summary_rounded_main_clean$`Performance metric` <-
  as.factor(
    c(
      "AUCTEST",
      "AUCDIFF",
      "TSSTEST",
      "TSSDIFF",
      "ORMIN-TEST",
      "ORMAX(SE+SP)-TEST",
      "ORMAX(SE+SP)-DIFF"
    )
  )
outerloop_cv_output_stats_summary_rounded_main_clean

write.csv(
  outerloop_cv_output_stats_summary_rounded_main_clean,
  file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/outerloop_cv_output_stats_summary_(estimated_generalization_error)_rounded_main_clean.csv",
  row.names = FALSE
)

# save all objects in R environment
save.image(
  "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil.RData"
)


# Run inner cross-validation procedure on full data set to tune final model ----

# create cross-validation folds (k = 4) for full data set
set.seed(25)
final_cv_folds <-
  get.checkerboard2(
    occs = data@coords[data@pa == 1,],
    envs = bio_env_c,
    bg = data@coords[data@pa == 0,],
    aggregation.factor = 4
  )

# train model via cross-validation
set.seed(25)
final_cv_model_untuned <-
  train(data = data,
        method = "Maxent",
        folds = final_cv_folds)

# plot the checkerboard partitions (k = 4) of the presence points
(
  final_cv_folds_presence_plot <-
    evalplot.grps(
      pts = data@coords[data@pa == 1,],
      pts.grp = final_cv_folds$occs.grp,
      envs = elevation_c,
      pts.size = 0.5
    )
)

# save the image output
ggsave(
  file = "final_cv_folds_presence_plot.png",
  final_cv_folds_presence_plot,
  device = "png",
  path = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil",
  width = 7,
  height = 7,
  units = "in",
  dpi = 400
)
ggsave(
  file = "final_cv_folds_presence_plot.pdf",
  final_cv_folds_presence_plot,
  device = "pdf",
  path = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil",
  width = 7,
  height = 7,
  units = "in"
)

# plot the checkerboard partitions (k = 4) of the background points
(
  final_cv_folds_background_plot <-
    evalplot.grps(
      pts = data@coords[data@pa == 0,],
      pts.grp = final_cv_folds$bg.grp,
      envs = elevation_c,
      pts.size = 0.5
    )
)

# save the image output
ggsave(
  file = "final_cv_folds_background_plot.png",
  final_cv_folds_background_plot,
  device = "png",
  path = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil",
  width = 7,
  height = 7,
  units = "in",
  dpi = 400
)
ggsave(
  file = "final_cv_folds_background_plot.pdf",
  final_cv_folds_background_plot,
  device = "pdf",
  path = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil",
  width = 7,
  height = 7,
  units = "in"
)

# save all objects in R environment
save.image(
  "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil.RData"
)


## Tune model hyperparameters ----

# see which hyperparameters can be tuned in the given model
getTunableArgs(final_cv_model_untuned)

# create a list of hyperparameter combinations to be tested, including L1 regularization multiplier values, feature classes, and number of iterations
hyperparam_combos <-
  list(
    reg = seq(0.1, 10, 0.1),
    fc = c("l", "lq", "lh", "lp", "lqp", "lqph"),
    iter = seq(500, 2000, 100)
  )
summary(hyperparam_combos)


### Optimize search (genetic algorithm) ----

# perform an optimize search to find the best possible combination of hyperparameters for the model chosen via a data-driven genetic algorithm from an initial random sample

{
  final_cv_optimizesearch_output_start_time <-
    proc.time()

  final_cv_optimizesearch_output <-
    optimizeModel(
      final_cv_model_untuned,
      hypers = hyperparam_combos,
      metric = "auc",
      pop = 30,
      gen = 5,
      interactive = FALSE,
      seed = 25
    )

  final_cv_optimizesearch_output_stop_time <-
    proc.time()
  final_cv_optimizesearch_output_run_time <-
    hms::hms(
      final_cv_optimizesearch_output_stop_time[["elapsed"]] - final_cv_optimizesearch_output_start_time[["elapsed"]]
    )
  print(str_c("Run time (HMS): ",
              final_cv_optimizesearch_output_run_time))
}

final_cv_optimizesearch_output
summary(final_cv_optimizesearch_output)

# save object to disk
save(final_cv_optimizesearch_output, file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/final_cv_optimizesearch_output.RData")

# check the output (ordered by test AUC in descending order by default), including only unique results
(
  final_cv_optimizesearch_output_results <-
    final_cv_optimizesearch_output@results %>% unique()
)

# plot the hyperparameter tuning search chart
(
  final_cv_optimizesearch_output_plot <-
    plot(final_cv_optimizesearch_output,
         title = "Hyperparameter Tuning:\nOptimize Search")
)

# save the image output
ggsave(
  file = "final_cv_optimizesearch_output_plot.pdf",
  final_cv_optimizesearch_output_plot,
  device = "pdf",
  path = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil",
  width = 7,
  height = 6,
  units = "in"
)
ggsave(
  file = str_c("final_cv_optimizesearch_output_plot.png"),
  final_cv_optimizesearch_output_plot,
  device = "png",
  path = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil",
  width = 7,
  height = 6,
  units = "in",
  dpi = 400
)

# save the best hyperparameter combination to an object
(final_hyperparam_optimizesearch <-
    final_cv_optimizesearch_output_results[1, ])

# save the best hyperparameter combination as a csv file
write.csv(final_hyperparam_optimizesearch,
          file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/final_hyperparam_optimizesearch.csv",
          row.names = FALSE)

# pull the best cv model from the tuning results and assign it to new model object
(final_cv_model_tuned <-
    final_cv_optimizesearch_output@models[[1]])

# save the output from the tuned final cv model as a txt file
final_cv_model_tuned
sink(
  "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/final_cv_model_tuned.txt"
)
print(final_cv_model_tuned)
sink()

final_cv_model_tuned_data <-
  final_cv_model_tuned@data

# save data from the tuned final cv model as a csv file
swd2csv(final_cv_model_tuned_data,
        file_name = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/final_cv_model_tuned_data.csv")

# save all objects in R environment
save.image(
  "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil.RData"
)


# Build and evaluate final model using full data set with tuned hyperparameters ----

# save the final optimized hyperparameter combination to an object
(final_hyperparam <-
   final_hyperparam_optimizesearch)

# save the final optimized hyperparameters as a csv file
write.csv(final_hyperparam,
          file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/final_hyperparam.csv",
          row.names = FALSE)

# build the final model using the full data set and with the optimal hyperparameter combo
set.seed(25)
final_model <-
  train(
    "Maxent",
    data = data,
    fc = final_hyperparam[, "fc"],
    reg = final_hyperparam[, "reg"],
    iter = final_hyperparam[, "iter"]
  )

# save the output of the final model as a txt file
final_model
sink(
  "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/final_model_output.txt"
)
print(final_model)
sink()

# save data from the final model as a csv file
swd2csv(final_model@data, file_name = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/final_model_data.csv")

# print the threshold values generated by maxent.jar
(ths_maxent_jar <- maxentTh(final_model))

# save output as a csv file
write.csv(ths_maxent_jar, file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/thresholds_values_maxent_jar.csv")

# compute the threshold values representing minimum training presence, equal training sensitivity and specificity, and maximum training sensitivity plus specificity, including the fractional predicted area and omission rate
(ths <- thresholds(final_model, type = "cloglog"))

# save output as a csv file
write.csv(ths, file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/thresholds_values.csv", row.names = FALSE)

# compute confusion matrix from final model
final_model_confusion_matrix <-
  confMatrix(final_model, type = "cloglog")

# save output as a csv file
write.csv(final_model_confusion_matrix,
          file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/final_model_confusion_matrix.csv",
          row.names = FALSE)

# create empty data frame for final model output
final_model_output_stats <-
  data.frame(matrix(NA, nrow = 1, ncol = 6))

colnames(final_model_output_stats) <-
  c("fc",
    "reg",
    "iter",
    "AUC",
    "TSS",
    "OR_max_train_se_plus_sp")
final_model_output_stats

# get optimal hyperparameter setting for feature classes
final_model_output_stats$fc <-
  final_hyperparam[, "fc"]

# get optimal hyperparameter setting for regularization multiplier
final_model_output_stats$reg <-
  final_hyperparam[, "reg"]

# get optimal hyperparameter setting for maximum number of iterations
final_model_output_stats$iter <-
  final_hyperparam[, "iter"]

# calculate train AUC
final_model_output_stats$AUC <-
  auc(final_model)

# calculate train TSS
final_model_output_stats$TSS <-
  tss(final_model)

# calculate train omission rate at maximum training sensitivity plus specificity threshold
final_model_output_stats$OR_max_train_se_plus_sp <-
  thresholds(final_model,
             type = "cloglog")[3, 4]

# check the final model output stats (AUC, TSS, and maxSSS threshold's OR)
(final_model_output_stats <- as_tibble(final_model_output_stats))

# compare with the outer loop output and summary stats
outerloop_cv_output_stats
outerloop_cv_output_stats_summary

str_c("Final model training AUC: ",
      final_model_output_stats$AUC %>%
        signif(4))
str_c(
  "[Estimated generalization error] Test AUC (mean +/- standard deviation): ",
  outerloop_cv_output_stats_summary %>%
    filter(performance_metric == "AUC") %>%
    dplyr::select(mean),
  " +/- ",
  outerloop_cv_output_stats_summary %>%
    filter(performance_metric == "AUC") %>%
    dplyr::select(sd)
)

str_c("Final model training TSS: ",
      final_model_output_stats$TSS %>%
        signif(4))
str_c(
  "[Estimated generalization error] Test TSS (mean +/- standard deviation): ",
  outerloop_cv_output_stats_summary %>%
    filter(performance_metric == "TSS") %>%
    dplyr::select(mean),
  " +/- ",
  outerloop_cv_output_stats_summary %>%
    filter(performance_metric == "TSS") %>%
    dplyr::select(sd)
)

str_c(
  "Final model training maxSSS OR: ",
  final_model_output_stats$OR_max_train_se_plus_sp %>%
    signif(4)
)
str_c(
  "[Estimated generalization error] Test maxSSS OR (mean +/- standard deviation): ",
  outerloop_cv_output_stats_summary %>%
    filter(performance_metric == "OR_max_train_se_plus_sp") %>%
    dplyr::select(mean),
  " +/- ",
  outerloop_cv_output_stats_summary %>%
    filter(performance_metric == "OR_max_train_se_plus_sp") %>%
    dplyr::select(sd)
)

# save the final model output stats as a csv file
write.csv(final_model_output_stats,
          file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/final_model_output_stats.csv",
          row.names = FALSE)

# make version of table with values rounded to 3 digits
final_model_output_stats_rounded <-
  final_model_output_stats %>%
  mutate(across(where(is.character), ~ as.factor(.x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x >= 0.001, round(.x, 3), .x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x < 0.001 & .x != 0, "<0.001", .x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x == 0, "0", .x)))
final_model_output_stats_rounded

write.csv(final_model_output_stats_rounded,
          file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/final_model_output_stats_rounded.csv",
          row.names = FALSE)

# save all objects in R environment
save.image(
  "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil.RData"
)


# Explore lambda parameters ----

# the lambda parameters (fitted lambda/weight coefficient value, minimum value, and maximum value) for each feature, including all raw variables and any versions of the variables transformed by different feature classes

# for final model
final_model_lambda_parameters <-
  separate(
    data = as.data.frame(final_model@model@lambdas),
    col = 1,
    into = c("feature", "lambda", "min", "max"),
    sep = ","
  ) %>% na.omit()

# for tuned outer loop models
outerloop_model_tuned_lambda_parameters_list <-
  sapply(1:num_outer_k_folds, function(x)
    NULL)

for (i in 1:num_outer_k_folds) {
  outerloop_model_tuned_lambda_parameters_list[[i]] <-
    separate(
      data = as.data.frame(outerloop_model_tuned_list[[i]]@model@lambdas),
      col = 1,
      into = c("feature", "lambda", "min", "max"),
      sep = ","
    ) %>% na.omit()
}

# quadratic (squared data values): [variable]^2
# product (product of two different variables): [feature 1]*[variable 2]
# forward hinge (rising linear slope function starting at the hinge value and clamped to 0 when less than the hinge value): '[variable]
# reverse hinge (negative linear slope, opposite of forward hinge): `[variable]
# threshold feature (0 when variable is less than the threshold and 1 otherwise: ([variable value]<[variable])

# non-zero coefficients (includes the following columns: feature used in the final model, including any variable transformations as indicated by "I()"; lambda coefficient value; minimum value, and maximum value)

# for final model
final_model_nonzero_coeff <- final_model@model@coeff

# for tuned outer loop models
outerloop_model_tuned_nonzero_coeff_list <-
  sapply(1:num_outer_k_folds, function(x)
    NULL)

for (i in 1:num_outer_k_folds) {
  outerloop_model_tuned_nonzero_coeff_list[[i]] <-
    outerloop_model_tuned_list[[i]]@model@coeff
}

# get number of non-zero coefficients (including those from feature transformations)

# for final model
final_model_nonzero_coeff_count <-
  as.data.frame(nrow(final_model@model@coeff))
colnames(final_model_nonzero_coeff_count) <-
  "nonzero_coeff_count"

# for tuned outer loop models
outerloop_model_tuned_nonzero_coeff_count_list <-
  sapply(1:num_outer_k_folds, function(x)
    NULL)

for (i in 1:num_outer_k_folds) {
  outerloop_model_tuned_nonzero_coeff_count_list[[i]] <-
    as.data.frame(nrow(outerloop_model_tuned_list[[i]]@model@coeff))
  colnames(outerloop_model_tuned_nonzero_coeff_count_list[[i]]) <-
    "nonzero_coeff_count"
}

# save the outputs as csv files

# for final model
write.csv(final_model_lambda_parameters,
          file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/final_model_lambda_parameters.csv",
          row.names = FALSE)
write.csv(final_model_nonzero_coeff,
          file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/final_model_nonzero_coeff.csv",
          row.names = FALSE)
write.csv(final_model_nonzero_coeff_count,
          file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/final_model_nonzero_coeff_count.csv",
          row.names = FALSE)

# for tuned outer loop models
for (i in 1:num_outer_k_folds) {
  write.csv(
    outerloop_model_tuned_lambda_parameters_list[[i]],
    file = str_c(
      "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/outerloop_model_tuned_lambda_parameters_",
      i,
      ".csv"
    ),
    row.names = FALSE
  )
  write.csv(
    outerloop_model_tuned_nonzero_coeff_list[[i]],
    file = str_c(
      "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/outerloop_model_tuned_nonzero_coeff_",
      i,
      ".csv"
    ),
    row.names = FALSE
  )
  write.csv(
    outerloop_model_tuned_nonzero_coeff_count_list[[i]],
    file = str_c(
      "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/outerloop_model_tuned_nonzero_coeff_count_",
      i,
      ".csv"
    ),
    row.names = FALSE
  )
}


# Plot response curves ----

if (dir.exists(
  "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/response_curves"
) == FALSE) {
  dir.create(
    "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/response_curves"
  )
}

if (dir.exists(
  "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/response_curves_rescaled"
) == FALSE) {
  dir.create(
    "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/response_curves_rescaled"
  )
}

if (dir.exists(
  "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/response_curves_plot_data"
) == FALSE) {
  dir.create(
    "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/response_curves_plot_data"
  )
}

# convert the final model and tuned outer loop models (SDMmodel objects) into dismo MaxEnt objects
final_model_dismo <- SDMmodel2MaxEnt(final_model)

outerloop_model_tuned_dismo_list <-
  sapply(1:num_outer_k_folds, function(x)
    NULL)

for (i in 1:num_outer_k_folds) {
  outerloop_model_tuned_dismo_list[[i]] <-
    SDMmodel2MaxEnt(outerloop_model_tuned_list[[i]])
}

# plot the marginal response curves for each variable (within the range of presences only) from the final model and tuned outer loop models; a rug of deciles is plotted on the horizonal axes
pdf(
  "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/response_curves/final_model_dismo_rc_margin_plot.pdf",
  width = 35,
  height = 35
)
response(final_model_dismo)
dev.off()

for (i in 1:num_outer_k_folds) {
  pdf(
    str_c(
      "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/response_curves/outerloop_model_tuned_dismo_rc_margin_plot_",
      i,
      ".pdf"
    ),
    width = 35,
    height = 35
  )
  response(outerloop_model_tuned_dismo_list[[i]])
  dev.off()
}

# plot the marginal response curves for each variable (within the range of presences only) from the final model and tuned outer loop models, this time via SDMtune; a rug of deciles is plotted on the horizonal axes, and the values are in cloglog scale; also, export the plot data (including the presence and background rugs) as csv files

{
  plan(multisession, workers = num_parallel_workers)
}  # initiate parallel processing

foreach(
  y = seq_along(bio_env_layer_names$variable_abbreviation),
  .options.future = list(seed = TRUE)
) %dofuture% {
  response_curve_marginal_plot <- plotResponse(
    final_model,
    var = bio_env_layer_names$variable_abbreviation[[y]],
    type = "cloglog",
    only_presence = TRUE,
    marginal = TRUE,
    rug = TRUE,
    color = "darkorange2"
  )
  ggsave(
    filename = str_c(
      "final_model_rc_margin_plot_",
      bio_env_layer_names$variable_abbreviation[[y]],
      ".png"
    ),
    plot = response_curve_marginal_plot,
    device = "png",
    path = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/response_curves/",
    width = 7,
    height = 6,
    units = "in",
    dpi = 400
  )
  ggsave(
    filename = str_c(
      "final_model_rc_margin_plot_",
      bio_env_layer_names$variable_abbreviation[[y]],
      ".pdf"
    ),
    plot = response_curve_marginal_plot,
    device = "pdf",
    path = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/response_curves/",
    width = 7,
    height = 6,
    units = "in"
  )

  response_curve_marginal_plot_rescaled <-
    response_curve_marginal_plot +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::xlab(if (bio_env_layer_names$units[[y]] != "") {
      str_c(bio_env_layer_names$variable_name[[y]],
            " (",
            bio_env_layer_names$units[[y]],
            ")")
    } else {
      str_c(bio_env_layer_names$variable_name[[y]])
    }) +
    ggplot2::ylab("Cloglog output (probability of presence)")
  ggsave(
    filename = str_c(
      "final_model_rc_margin_plot_",
      bio_env_layer_names$variable_abbreviation[[y]],
      ".png"
    ),
    plot = response_curve_marginal_plot_rescaled,
    device = "png",
    path = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/response_curves_rescaled/",
    width = 7,
    height = 6,
    units = "in",
    dpi = 400
  )
  ggsave(
    filename = str_c(
      "final_model_rc_margin_plot_",
      bio_env_layer_names$variable_abbreviation[[y]],
      ".pdf"
    ),
    plot = response_curve_marginal_plot_rescaled,
    device = "pdf",
    path = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/response_curves_rescaled/",
    width = 7,
    height = 6,
    units = "in"
  )

  response_curve_marginal_plot_rescaled_data <-
    ggplot2::ggplot_build(response_curve_marginal_plot_rescaled)$data[[1]]
  response_curve_marginal_plot_rescaled_data_rugP <-
    ggplot2::ggplot_build(response_curve_marginal_plot_rescaled)$data[[2]]
  response_curve_marginal_plot_rescaled_data_rugBG <-
    ggplot2::ggplot_build(response_curve_marginal_plot_rescaled)$data[[3]]
  write.csv(
    response_curve_marginal_plot_rescaled_data,
    file = str_c(
      "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/response_curves_plot_data/",
      "final_model_rc_margin_plot_",
      bio_env_layer_names$variable_abbreviation[[y]],
      "_data",
      ".csv"
    ),
    row.names = FALSE
  )
  write.csv(
    response_curve_marginal_plot_rescaled_data_rugP,
    file = str_c(
      "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/response_curves_plot_data/",
      "final_model_rc_margin_plot_",
      bio_env_layer_names$variable_abbreviation[[y]],
      "_rugP",
      ".csv"
    ),
    row.names = FALSE
  )
  write.csv(
    response_curve_marginal_plot_rescaled_data_rugBG,
    file = str_c(
      "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/response_curves_plot_data/",
      "final_model_rc_margin_plot_",
      bio_env_layer_names$variable_abbreviation[[y]],
      "_rugBG",
      ".csv"
    ),
    row.names = FALSE
  )
}

foreach(i = 1:num_outer_k_folds,
        .options.future = list(seed = TRUE)) %dofuture% {
          foreach(
            y = seq_along(bio_env_layer_names$variable_abbreviation),
            .options.future = list(seed = TRUE)
          ) %dofuture% {
            response_curve_marginal_plot <- plotResponse(
              outerloop_model_tuned_list[[i]],
              var = bio_env_layer_names$variable_abbreviation[[y]],
              type = "cloglog",
              only_presence = TRUE,
              marginal = TRUE,
              rug = TRUE,
              color = "darkorange2"
            )
            ggsave(
              filename = str_c(
                "outerloop_model_tuned_rc_margin_plot_",
                bio_env_layer_names$variable_abbreviation[[y]],
                "_",
                i,
                ".png"
              ),
              plot = response_curve_marginal_plot,
              device = "png",
              path = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/response_curves/",
              width = 7,
              height = 6,
              units = "in",
              dpi = 400
            )
            ggsave(
              filename = str_c(
                "outerloop_model_tuned_rc_margin_plot_",
                bio_env_layer_names$variable_abbreviation[[y]],
                "_",
                i,
                ".pdf"
              ),
              plot = response_curve_marginal_plot,
              device = "pdf",
              path = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/response_curves/",
              width = 7,
              height = 6,
              units = "in"
            )

            response_curve_marginal_plot_rescaled <-
              response_curve_marginal_plot +
              ggplot2::scale_y_continuous(limits = c(0, 1)) +
              ggplot2::xlab(if (bio_env_layer_names$units[[y]] != "") {
                str_c(bio_env_layer_names$variable_name[[y]],
                      " (",
                      bio_env_layer_names$units[[y]],
                      ")")
              } else {
                str_c(bio_env_layer_names$variable_name[[y]])
              }) +
              ggplot2::ylab("Cloglog output (probability of presence)")
            ggsave(
              filename = str_c(
                "outerloop_model_rc_margin_plot_",
                bio_env_layer_names$variable_abbreviation[[y]],
                "_",
                i,
                ".png"
              ),
              plot = response_curve_marginal_plot_rescaled,
              device = "png",
              path = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/response_curves_rescaled/",
              width = 7,
              height = 6,
              units = "in",
              dpi = 400
            )
            ggsave(
              filename = str_c(
                "outerloop_model_rc_margin_plot_",
                bio_env_layer_names$variable_abbreviation[[y]],
                "_",
                i,
                ".pdf"
              ),
              plot = response_curve_marginal_plot_rescaled,
              device = "pdf",
              path = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/response_curves_rescaled/",
              width = 7,
              height = 6,
              units = "in"
            )

            response_curve_marginal_plot_rescaled_data <-
              ggplot2::ggplot_build(response_curve_marginal_plot_rescaled)$data[[1]]
            response_curve_marginal_plot_rescaled_data_rugP <-
              ggplot2::ggplot_build(response_curve_marginal_plot_rescaled)$data[[2]]
            response_curve_marginal_plot_rescaled_data_rugBG <-
              ggplot2::ggplot_build(response_curve_marginal_plot_rescaled)$data[[3]]
            write.csv(
              response_curve_marginal_plot_rescaled_data,
              file = str_c(
                "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/response_curves_plot_data/",
                "outerloop_model_tuned_rc_margin_plot_",
                bio_env_layer_names$variable_abbreviation[[y]],
                "_",
                i,
                "_data",
                ".csv"
              ),
              row.names = FALSE
            )
            write.csv(
              response_curve_marginal_plot_rescaled_data_rugP,
              file = str_c(
                "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/response_curves_plot_data/",
                "outerloop_model_tuned_rc_margin_plot_",
                bio_env_layer_names$variable_abbreviation[[y]],
                "_",
                i,
                "_rugP",
                ".csv"
              ),
              row.names = FALSE
            )
            write.csv(
              response_curve_marginal_plot_rescaled_data_rugBG,
              file = str_c(
                "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/response_curves_plot_data/",
                "outerloop_model_tuned_rc_margin_plot_",
                bio_env_layer_names$variable_abbreviation[[y]],
                "_",
                i,
                "_rugBG",
                ".csv"
              ),
              row.names = FALSE
            )
          }
        }

{
  plan(sequential)
}  # end parallel processing

# plot the univariate response curves for each variable (within the range of presences only) from the final model and tuned outer loop models, this time via SDMtune; a rug of deciles is plotted on the horizonal axes, and the values are in cloglog scale; also, export the plot data (including the presence and background rugs) as csv files

{
  plan(multisession, workers = num_parallel_workers)
}  # initiate parallel processing

foreach(
  y = seq_along(bio_env_layer_names$variable_abbreviation),
  .options.future = list(seed = TRUE)
) %dofuture% {
  response_curve_univariate_plot <- plotResponse(
    final_model,
    var = bio_env_layer_names$variable_abbreviation[[y]],
    type = "cloglog",
    only_presence = TRUE,
    marginal = FALSE,
    rug = TRUE,
    color = "blue"
  )
  ggsave(
    filename = str_c(
      "final_model_rc_univar_plot_",
      bio_env_layer_names$variable_abbreviation[[y]],
      ".png"
    ),
    plot = response_curve_univariate_plot,
    device = "png",
    path = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/response_curves/",
    width = 7,
    height = 6,
    units = "in",
    dpi = 400
  )
  ggsave(
    filename = str_c(
      "final_model_rc_univar_plot_",
      bio_env_layer_names$variable_abbreviation[[y]],
      ".pdf"
    ),
    plot = response_curve_univariate_plot,
    device = "pdf",
    path = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/response_curves/",
    width = 7,
    height = 6,
    units = "in"
  )

  response_curve_univariate_plot_rescaled <-
    response_curve_univariate_plot +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::xlab(if (bio_env_layer_names$units[[y]] != "") {
      str_c(bio_env_layer_names$variable_name[[y]],
            " (",
            bio_env_layer_names$units[[y]],
            ")")
    } else {
      str_c(bio_env_layer_names$variable_name[[y]])
    }) +
    ggplot2::ylab("Cloglog output (probability of presence)")
  ggsave(
    filename = str_c(
      "final_model_rc_univar_plot_",
      bio_env_layer_names$variable_abbreviation[[y]],
      ".png"
    ),
    plot = response_curve_univariate_plot_rescaled,
    device = "png",
    path = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/response_curves_rescaled/",
    width = 7,
    height = 6,
    units = "in",
    dpi = 400
  )
  ggsave(
    filename = str_c(
      "final_model_rc_univar_plot_",
      bio_env_layer_names$variable_abbreviation[[y]],
      ".pdf"
    ),
    plot = response_curve_univariate_plot_rescaled,
    device = "pdf",
    path = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/response_curves_rescaled/",
    width = 7,
    height = 6,
    units = "in"
  )

  response_curve_univariate_plot_rescaled_data <-
    ggplot2::ggplot_build(response_curve_univariate_plot_rescaled)$data[[1]]
  response_curve_univariate_plot_rescaled_data_rugP <-
    ggplot2::ggplot_build(response_curve_univariate_plot_rescaled)$data[[2]]
  response_curve_univariate_plot_rescaled_data_rugBG <-
    ggplot2::ggplot_build(response_curve_univariate_plot_rescaled)$data[[3]]
  write.csv(
    response_curve_univariate_plot_rescaled_data,
    file = str_c(
      "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/response_curves_plot_data/",
      "final_model_rc_univar_plot_",
      bio_env_layer_names$variable_abbreviation[[y]],
      "_data",
      ".csv"
    ),
    row.names = FALSE
  )
  write.csv(
    response_curve_univariate_plot_rescaled_data_rugP,
    file = str_c(
      "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/response_curves_plot_data/",
      "final_model_rc_univar_plot_",
      bio_env_layer_names$variable_abbreviation[[y]],
      "_rugP",
      ".csv"
    ),
    row.names = FALSE
  )
  write.csv(
    response_curve_univariate_plot_rescaled_data_rugBG,
    file = str_c(
      "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/response_curves_plot_data/",
      "final_model_rc_univar_plot_",
      bio_env_layer_names$variable_abbreviation[[y]],
      "_rugBG",
      ".csv"
    ),
    row.names = FALSE
  )
}

foreach(i = 1:num_outer_k_folds,
        .options.future = list(seed = TRUE)) %dofuture% {
          foreach(
            y = seq_along(bio_env_layer_names$variable_abbreviation),
            .options.future = list(seed = TRUE)
          ) %dofuture% {
            response_curve_univariate_plot <- plotResponse(
              outerloop_model_tuned_list[[i]],
              var = bio_env_layer_names$variable_abbreviation[[y]],
              type = "cloglog",
              only_presence = TRUE,
              marginal = FALSE,
              rug = TRUE,
              color = "blue"
            )
            ggsave(
              filename = str_c(
                "outerloop_model_tuned_rc_univar_plot_",
                bio_env_layer_names$variable_abbreviation[[y]],
                "_",
                i,
                ".png"
              ),
              plot = response_curve_univariate_plot,
              device = "png",
              path = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/response_curves/",
              width = 7,
              height = 6,
              units = "in",
              dpi = 400
            )
            ggsave(
              filename = str_c(
                "outerloop_model_tuned_rc_univar_plot_",
                bio_env_layer_names$variable_abbreviation[[y]],
                "_",
                i,
                ".pdf"
              ),
              plot = response_curve_univariate_plot,
              device = "pdf",
              path = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/response_curves/",
              width = 7,
              height = 6,
              units = "in"
            )

            response_curve_univariate_plot_rescaled <-
              response_curve_univariate_plot +
              ggplot2::scale_y_continuous(limits = c(0, 1)) +
              ggplot2::xlab(if (bio_env_layer_names$units[[y]] != "") {
                str_c(bio_env_layer_names$variable_name[[y]],
                      " (",
                      bio_env_layer_names$units[[y]],
                      ")")
              } else {
                str_c(bio_env_layer_names$variable_name[[y]])
              }) +
              ggplot2::ylab("Cloglog output (probability of presence)")
            ggsave(
              filename = str_c(
                "outerloop_model_rc_univar_plot_",
                bio_env_layer_names$variable_abbreviation[[y]],
                "_",
                i,
                ".png"
              ),
              plot = response_curve_univariate_plot_rescaled,
              device = "png",
              path = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/response_curves_rescaled/",
              width = 7,
              height = 6,
              units = "in",
              dpi = 400
            )
            ggsave(
              filename = str_c(
                "outerloop_model_rc_univar_plot_",
                bio_env_layer_names$variable_abbreviation[[y]],
                "_",
                i,
                ".pdf"
              ),
              plot = response_curve_univariate_plot_rescaled,
              device = "pdf",
              path = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/response_curves_rescaled/",
              width = 7,
              height = 6,
              units = "in"
            )

            response_curve_univariate_plot_rescaled_data <-
              ggplot2::ggplot_build(response_curve_univariate_plot_rescaled)$data[[1]]
            response_curve_univariate_plot_rescaled_data_rugP <-
              ggplot2::ggplot_build(response_curve_univariate_plot_rescaled)$data[[2]]
            response_curve_univariate_plot_rescaled_data_rugBG <-
              ggplot2::ggplot_build(response_curve_univariate_plot_rescaled)$data[[3]]
            write.csv(
              response_curve_univariate_plot_rescaled_data,
              file = str_c(
                "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/response_curves_plot_data/",
                "outerloop_model_tuned_rc_univar_plot_",
                bio_env_layer_names$variable_abbreviation[[y]],
                "_",
                i,
                "_data",
                ".csv"
              ),
              row.names = FALSE
            )
            write.csv(
              response_curve_univariate_plot_rescaled_data_rugP,
              file = str_c(
                "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/response_curves_plot_data/",
                "outerloop_model_tuned_rc_univar_plot_",
                bio_env_layer_names$variable_abbreviation[[y]],
                "_",
                i,
                "_rugP",
                ".csv"
              ),
              row.names = FALSE
            )
            write.csv(
              response_curve_univariate_plot_rescaled_data_rugBG,
              file = str_c(
                "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/response_curves_plot_data/",
                "outerloop_model_tuned_rc_univar_plot_",
                bio_env_layer_names$variable_abbreviation[[y]],
                "_",
                i,
                "_rugBG",
                ".csv"
              ),
              row.names = FALSE
            )
          }
        }

{
  plan(sequential)
}  # end parallel processing

# save all objects in R environment
save.image(
  "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil.RData"
)


# Assess multicollinearity and relative importance of environmental variables used in final models ----

## Permutation importance ----

# check Maxent-generated relative permutation importance (variable contribution, or impact each variable had on predicting distribution) for outer loop tuned models (with values from each model averaged together) as well as final models (tuned cross-validation model and full non-cross-validation model), ordered by permutation importance

### Outer loop tuned models (averaged) ----

# get relative variable contribution rankings as calculated by Maxent, sorting results by permutation importance in descending order; exclude percent contribution values (which, unlike permutation importance, are affected by the order in which Maxent evaluated each variable during the modeling process) from the resulting data frame
outerloop_model_tuned_var_imp_list <-
  sapply(1:num_outer_k_folds, function(x)
    NULL)

for (i in 1:num_outer_k_folds) {
  outerloop_model_tuned_var_imp_list[[i]] <-
    maxentVarImp(outerloop_model_tuned_list[[i]]) %>%
    arrange(desc(Permutation_importance)) %>%
    dplyr::select(-Percent_contribution)
}

# combine list's data frames via union join; save output as a csv file
outerloop_model_tuned_var_imp_joined_df <-
  outerloop_model_tuned_var_imp_list %>%
  reduce(union_all) %>%
  tibble() %>%
  arrange(desc(Permutation_importance)) %>%
  left_join(
    dplyr::select(bio_env_layer_names, variable_name, variable_abbreviation),
    by = join_by("Variable" == "variable_abbreviation")
  ) %>%
  relocate(Variable_description = variable_name, .after = Variable)
outerloop_model_tuned_var_imp_joined_df

write.csv(outerloop_model_tuned_var_imp_joined_df,
          file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/outerloop_model_tuned_var_imp_joined_df.csv",
          row.names = FALSE)

# calculate summary statistics for relative permutation importance from outer loop tuned models (where sample size (n) = number of outer loop cross-validation folds (k)); these values, particularly the mean and standard deviation, can then be reported alongside the final model as estimated generalization error (estimated average relative permutation importance of variables used in final model); save output as a csv file
outerloop_model_tuned_var_imp_summary <-
  outerloop_model_tuned_var_imp_list %>%
  reduce(union_all) %>%
  group_by(Variable) %>%
  summarize(tibble(
    across(where(is.numeric), ~ mean(.x, na.rm = TRUE), .names = "mean_{.col}"),
    across(where(is.numeric), ~ sd(.x, na.rm = TRUE), .names = "sd_{.col}"),
    across(where(is.numeric), ~ min(.x, na.rm = TRUE), .names = "min_{.col}"),
    across(where(is.numeric), ~ max(.x, na.rm = TRUE), .names = "max_{.col}"),
    across(
      where(is.numeric),
      ~ quantile(.x, probs = 0.25, na.rm = TRUE),
      .names = "q1_{.col}"
    ),
    # 1st quantile
    across(where(is.numeric), ~ median(.x, na.rm = TRUE), .names = "median_{.col}"),
    across(
      where(is.numeric),
      ~ quantile(.x, probs = 0.75, na.rm = TRUE),
      .names = "q3_{.col}"
    ),
    # 3rd quantile
    across(where(is.numeric), ~ IQR(.x, na.rm = TRUE), .names = "IQR_{.col}")
  )) %>%
  arrange(desc(mean_Permutation_importance)) %>%
  left_join(
    dplyr::select(bio_env_layer_names, variable_name, variable_abbreviation),
    by = join_by("Variable" == "variable_abbreviation")
  ) %>%
  relocate(Variable_description = variable_name, .after = Variable)
outerloop_model_tuned_var_imp_summary

write.csv(outerloop_model_tuned_var_imp_summary,
          file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/outerloop_model_tuned_var_imp_summary.csv",
          row.names = FALSE)

# make version of table with values rounded to 3 digits
outerloop_model_tuned_var_imp_summary_rounded <-
  outerloop_model_tuned_var_imp_summary %>%
  mutate(across(where(is.character), ~ as.factor(.x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x >= 0.001, round(.x, 3), .x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x < 0.001 & .x != 0, "<0.001", .x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x == 0, "0", .x)))

write.csv(
  outerloop_model_tuned_var_imp_summary_rounded,
  file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/outerloop_model_tuned_var_imp_summary_rounded.csv",
  row.names = FALSE
)

# show same data frame, but this time with just the most relevant/main statistics (mean, standard deviation, and min/max), then save output as a csv file
outerloop_model_tuned_var_imp_main_stats <-
  outerloop_model_tuned_var_imp_summary %>% dplyr::select(
    Variable,
    Variable_description,
    mean_Permutation_importance,
    sd_Permutation_importance,
    min_Permutation_importance,
    max_Permutation_importance
  )
outerloop_model_tuned_var_imp_main_stats

write.csv(outerloop_model_tuned_var_imp_main_stats,
          file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/outerloop_model_tuned_var_imp_main_stats.csv",
          row.names = FALSE)

# make version of table with values rounded to 3 digits
outerloop_model_tuned_var_imp_main_stats_rounded <-
  outerloop_model_tuned_var_imp_main_stats %>%
  mutate(across(where(is.character), ~ as.factor(.x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x >= 0.001, round(.x, 3), .x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x < 0.001 & .x != 0, "<0.001", .x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x == 0, "0", .x)))

write.csv(
  outerloop_model_tuned_var_imp_main_stats_rounded,
  file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/outerloop_model_tuned_var_imp_main_stats_rounded.csv",
  row.names = FALSE
)

# adjust/transform summarized and individual data to prepare for plotting
outerloop_model_tuned_var_imp_main_stats_plot_data <-
  outerloop_model_tuned_var_imp_main_stats %>%
  arrange(mean_Permutation_importance)
outerloop_model_tuned_var_imp_joined_df_plot_data <-
  outerloop_model_tuned_var_imp_joined_df %>%
  arrange(Permutation_importance)

outerloop_model_tuned_var_imp_main_stats_plot_data$Variable <-
  factor(outerloop_model_tuned_var_imp_main_stats_plot_data$Variable,
         levels = outerloop_model_tuned_var_imp_main_stats_plot_data$Variable)
outerloop_model_tuned_var_imp_joined_df_plot_data$Variable <-
  factor(
    outerloop_model_tuned_var_imp_joined_df_plot_data$Variable,
    levels = outerloop_model_tuned_var_imp_joined_df_plot_data$Variable %>% unique()
  )

outerloop_model_tuned_var_imp_main_stats_plot_data$Variable_description <-
  factor(
    outerloop_model_tuned_var_imp_main_stats_plot_data$Variable_description,
    levels = outerloop_model_tuned_var_imp_main_stats_plot_data$Variable_description
  )
outerloop_model_tuned_var_imp_joined_df_plot_data$Variable_description <-
  factor(
    outerloop_model_tuned_var_imp_joined_df_plot_data$Variable_description,
    levels = outerloop_model_tuned_var_imp_joined_df_plot_data$Variable_description %>% unique()
  )

outerloop_model_tuned_var_imp_main_stats_plot_data <-
  outerloop_model_tuned_var_imp_main_stats_plot_data %>%
  group_by(Variable_description) %>%
  summarize(across(where(is.numeric), ~ `/`(.x, 100)))
outerloop_model_tuned_var_imp_joined_df_plot_data <-
  outerloop_model_tuned_var_imp_joined_df_plot_data %>%
  group_by(Variable_description) %>%
  mutate(Permutation_importance = Permutation_importance / 100, .keep = "none")

# plot bar graph showing mean permutation importance with individual sample data points added for each variable, then save the image output
(
  outerloop_model_tuned_var_imp_bargraph <-
    ggplot(
      data = outerloop_model_tuned_var_imp_main_stats_plot_data,
      aes(x = Variable_description, y = mean_Permutation_importance)
    ) +
    labs(x = "Predictors",
         y = "Mean permutation importance") +
    scale_y_continuous(labels = scales::label_percent()) +
    coord_flip() +
    geom_bar(
      position = "dodge",
      stat = "identity",
      fill = alpha(custom_color_pal_1$green2_zomp, 0.4),
      color = alpha(custom_color_pal_1$gray_paynes, 0.9),
      linewidth = 0.6
    ) +
    geom_point(
      data = (
        outerloop_model_tuned_var_imp_joined_df_plot_data %>% filter(max(Permutation_importance) > 0)
      ),
      aes(x = Variable_description, y = Permutation_importance),
      size = 2.3,
      fill = alpha(custom_color_pal_1$yellow_saffron, 0.8),
      color = alpha(custom_color_pal_1$orange1_crayola, 0.9),
      shape = 23,
      stroke = 1.1
    ) +  # for variables with at least one permutation importance sample value greater than zero, add individual data points (n = 5)
    theme_bw() +
    theme(axis.title.x = element_text(
      color = "#444444",
      size = 24,
      margin = margin(20, 20, 20, 20)
    )) +
    theme(axis.title.y = element_text(
      color = "#444444",
      size = 24,
      margin = margin(20, 20, 20, 20)
    )) +
    theme(axis.text = element_text(
      color = "#444444", size = 16
    ))
)

ggsave(
  file = "outerloop_model_tuned_var_imp_bargraph.pdf",
  outerloop_model_tuned_var_imp_bargraph,
  device = "pdf",
  path = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil",
  width = 20,
  height = 28,
  units = "in"
)
ggsave(
  file = "outerloop_model_tuned_var_imp_bargraph.png",
  outerloop_model_tuned_var_imp_bargraph,
  device = "png",
  path = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil",
  width = 20,
  height = 28,
  units = "in",
  dpi = 400
)


### Final tuned model ----

# get relative variable contribution rankings as calculated by Maxent, sorting results by permutation importance in descending order; exclude percent contribution values (which, unlike permutation importance, are affected by the order in which Maxent evaluated each variable during the modeling process) from the resulting data frame
final_model_var_imp <-
  maxentVarImp(final_model) %>%
  tibble() %>%
  arrange(desc(Permutation_importance)) %>%
  dplyr::select(-Percent_contribution) %>%
  left_join(
    dplyr::select(bio_env_layer_names, variable_name, variable_abbreviation),
    by = join_by("Variable" == "variable_abbreviation")
  ) %>%
  relocate(Variable_description = variable_name, .after = Variable)
final_model_var_imp

# save output as a csv file
write.csv(final_model_var_imp,
          file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/final_model_var_imp.csv",
          row.names = FALSE)

# make version of table with values rounded to 3 digits
final_model_var_imp_rounded <-
  final_model_var_imp %>%
  mutate(across(where(is.character), ~ as.factor(.x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x >= 0.001, round(.x, 3), .x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x < 0.001 & .x != 0, "<0.001", .x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x == 0, "0", .x)))

write.csv(final_model_var_imp_rounded,
          file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/final_model_var_imp_rounded.csv",
          row.names = FALSE)

# adjust/transform summarized and individual data to prepare for plotting
final_model_var_imp_plot_data <-
  final_model_var_imp %>%
  arrange(Permutation_importance)

final_model_var_imp_plot_data$Variable <-
  factor(final_model_var_imp_plot_data$Variable,
         levels = final_model_var_imp_plot_data$Variable)
final_model_var_imp_plot_data$Variable_description <-
  factor(
    final_model_var_imp_plot_data$Variable_description,
    levels = final_model_var_imp_plot_data$Variable_description
  )

final_model_var_imp_plot_data <-
  final_model_var_imp_plot_data %>%
  group_by(Variable_description) %>%
  mutate(Permutation_importance = Permutation_importance / 100, .keep = "none")

# plot bar graph showing permutation importance in final model, then save the image output
(
  final_model_var_imp_bargraph <-
    ggplot(data = final_model_var_imp_plot_data,
           aes(x = Variable_description, y = Permutation_importance)) +
    labs(x = "Predictors", y = "Permutation importance") +
    scale_y_continuous(labels = scales::label_percent()) +
    coord_flip() +
    geom_bar(
      position = "dodge",
      stat = "identity",
      fill = alpha(custom_color_pal_1$green2_zomp, 0.4),
      color = alpha(custom_color_pal_1$gray_paynes, 0.9),
      linewidth = 0.6
    ) +
    theme_bw() +
    theme(axis.title.x = element_text(
      color = "#444444",
      size = 24,
      margin = margin(20, 20, 20, 20)
    )) +
    theme(axis.title.y = element_text(
      color = "#444444",
      size = 24,
      margin = margin(20, 20, 20, 20)
    )) +
    theme(axis.text = element_text(
      color = "#444444", size = 16
    ))
)

ggsave(
  file = "final_model_var_imp_bargraph.pdf",
  final_model_var_imp_bargraph,
  device = "pdf",
  path = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil",
  width = 20,
  height = 28,
  units = "in"
)
ggsave(
  file = "final_model_var_imp_bargraph.png",
  final_model_var_imp_bargraph,
  device = "png",
  path = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil",
  width = 20,
  height = 28,
  units = "in",
  dpi = 400
)

# save all objects in R environment
save.image(
  "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil.RData"
)


## Jackknife testing for variable importance ----

# check variable importance via jackknife testing (alternative method for calculating variable contribution, or impact each variable had on predicting distribution, by measuring effect on model when excluding one variable at a time as well as for each variable in isolation) for outer loop tuned models (with values from each model averaged together) as well as final models (tuned cross-validation model and full non-cross-validation model)

### Outer loop tuned models (averaged) ----

# get jackknife-tested variable importance from each outer loop tuned model, using corresponding outer loop withheld test data sets to calculate test AUC for model built without variable and then with variable in isolation
outerloop_model_tuned_var_imp_jk_list <-
  sapply(1:num_outer_k_folds, function(x)
    NULL)

{
  outerloop_model_tuned_var_imp_jk_list_start_time <-
    proc.time()

  {
    plan(multisession, workers = num_parallel_workers)
    }  # initiate parallel processing

  outerloop_model_tuned_var_imp_jk_list <-
    foreach(
      i = 1:num_outer_k_folds,
      .options.future = list(seed = TRUE,
                             packages = c("SDMtune"))
    ) %dofuture% {
      set.seed(25)
      outerloop_model_tuned_var_imp_jk <-
        doJk(
          outerloop_model_tuned_list[[i]],
          metric = "auc",
          test = outerloop_test_subset_list[[i]],
          with_only = TRUE
        )
      outerloop_model_tuned_var_imp_jk
    }

  write.csv(final_model_var_imp,
            file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/final_model_var_imp.csv",
            row.names = FALSE)

  {
    plan(sequential)
  }  # end parallel processing

  outerloop_model_tuned_var_imp_jk_list_stop_time <-
    proc.time()
  outerloop_model_tuned_var_imp_jk_list_run_time <-
    hms::hms(
      outerloop_model_tuned_var_imp_jk_list_stop_time[["elapsed"]] - outerloop_model_tuned_var_imp_jk_list_start_time[["elapsed"]]
    )
  print(str_c(
    "Run time (HMS): ",
    outerloop_model_tuned_var_imp_jk_list_run_time
  ))
  }

# combine list's data frames via union join; save output as a csv file
outerloop_model_tuned_var_imp_jk_full_joined_df <-
  outerloop_model_tuned_var_imp_jk_list %>%
  reduce(union_all) %>%
  tibble() %>%
  arrange(desc(Test_AUC_withonly)) %>%
  left_join(
    dplyr::select(bio_env_layer_names, variable_name, variable_abbreviation),
    by = join_by("Variable" == "variable_abbreviation")
  ) %>%
  relocate(Variable_description = variable_name, .after = Variable)
outerloop_model_tuned_var_imp_jk_full_joined_df

write.csv(
  outerloop_model_tuned_var_imp_jk_full_joined_df,
  file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/outerloop_model_tuned_var_imp_jk_full_joined_df.csv",
  row.names = FALSE
)

# make version of combined data frames using just the "withonly" test AUC metric; save output as a csv file
outerloop_model_tuned_var_imp_jk_test_AUC_withonly_joined_df <-
  outerloop_model_tuned_var_imp_jk_full_joined_df %>% dplyr::select(Variable, Variable_description, Test_AUC_withonly)
outerloop_model_tuned_var_imp_jk_test_AUC_withonly_joined_df

write.csv(
  outerloop_model_tuned_var_imp_jk_test_AUC_withonly_joined_df,
  file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/outerloop_model_tuned_var_imp_jk_test_AUC_withonly_joined_df.csv",
  row.names = FALSE
)

# calculate summary statistics for jackknife-tested variable importance from outer loop tuned models (where sample size (n) = number of outer loop cross-validation folds (k)); these values, particularly the mean and standard deviation, can then be reported alongside the final model as estimated generalization error (estimated average jackknife-tested variable importance in final model tested); sort by "withonly" test AUC metric, then save output as a csv file
outerloop_model_tuned_var_imp_jk_full_summary <-
  outerloop_model_tuned_var_imp_jk_list %>%
  reduce(union_all) %>%
  group_by(Variable) %>%
  summarize(tibble(
    across(where(is.numeric), ~ mean(.x, na.rm = TRUE), .names = "mean_{.col}"),
    across(where(is.numeric), ~ sd(.x, na.rm = TRUE), .names = "sd_{.col}"),
    across(where(is.numeric), ~ min(.x, na.rm = TRUE), .names = "min_{.col}"),
    across(where(is.numeric), ~ max(.x, na.rm = TRUE), .names = "max_{.col}"),
    across(
      where(is.numeric),
      ~ quantile(.x, probs = 0.25, na.rm = TRUE),
      .names = "q1_{.col}"
    ),
    # 1st quantile
    across(where(is.numeric), ~ median(.x, na.rm = TRUE), .names = "median_{.col}"),
    across(
      where(is.numeric),
      ~ quantile(.x, probs = 0.75, na.rm = TRUE),
      .names = "q3_{.col}"
    ),
    # 3rd quantile
    across(where(is.numeric), ~ IQR(.x, na.rm = TRUE), .names = "IQR_{.col}")
  )) %>%
  arrange(desc(mean_Test_AUC_withonly)) %>%
  left_join(
    dplyr::select(bio_env_layer_names, variable_name, variable_abbreviation),
    by = join_by("Variable" == "variable_abbreviation")
  ) %>%
  relocate(Variable_description = variable_name, .after = Variable)
outerloop_model_tuned_var_imp_jk_full_summary

write.csv(
  outerloop_model_tuned_var_imp_jk_full_summary,
  file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/outerloop_model_tuned_var_imp_jk_full_summary.csv",
  row.names = FALSE
)

# make version of table with values rounded to 3 digits
outerloop_model_tuned_var_imp_jk_full_summary_rounded <-
  outerloop_model_tuned_var_imp_jk_full_summary %>%
  mutate(across(where(is.character), ~ as.factor(.x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x >= 0.001, round(.x, 3), .x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x < 0.001 & .x != 0, "<0.001", .x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x == 0, "0", .x)))

# save output as a csv file
write.csv(
  outerloop_model_tuned_var_imp_jk_full_summary_rounded,
  file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/outerloop_model_tuned_var_imp_jk_full_summary_rounded.csv",
  row.names = FALSE
)

# make version using just the "withonly" test AUC metric
outerloop_model_tuned_var_imp_jk_test_AUC_withonly_summary <-
  outerloop_model_tuned_var_imp_jk_list %>%
  reduce(union_all) %>%
  dplyr::select(Variable, Test_AUC_withonly) %>%
  group_by(Variable) %>%
  summarize(tibble(
    across(where(is.numeric), ~ mean(.x, na.rm = TRUE), .names = "mean_{.col}"),
    across(where(is.numeric), ~ sd(.x, na.rm = TRUE), .names = "sd_{.col}"),
    across(where(is.numeric), ~ min(.x, na.rm = TRUE), .names = "min_{.col}"),
    across(where(is.numeric), ~ max(.x, na.rm = TRUE), .names = "max_{.col}"),
    across(
      where(is.numeric),
      ~ quantile(.x, probs = 0.25, na.rm = TRUE),
      .names = "q1_{.col}"
    ),
    # 1st quantile
    across(where(is.numeric), ~ median(.x, na.rm = TRUE), .names = "median_{.col}"),
    across(
      where(is.numeric),
      ~ quantile(.x, probs = 0.75, na.rm = TRUE),
      .names = "q3_{.col}"
    ),
    # 3rd quantile
    across(where(is.numeric), ~ IQR(.x, na.rm = TRUE), .names = "IQR_{.col}")
  )) %>%
  arrange(desc(mean_Test_AUC_withonly)) %>%
  left_join(
    dplyr::select(bio_env_layer_names, variable_name, variable_abbreviation),
    by = join_by("Variable" == "variable_abbreviation")
  ) %>%
  relocate(Variable_description = variable_name, .after = Variable)
outerloop_model_tuned_var_imp_jk_test_AUC_withonly_summary

# save output as a csv file
write.csv(
  outerloop_model_tuned_var_imp_jk_test_AUC_withonly_summary,
  file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/outerloop_model_tuned_var_imp_jk_test_AUC_withonly_summary.csv",
  row.names = FALSE
)

# make version of table with values rounded to 3 digits
outerloop_model_tuned_var_imp_jk_test_AUC_withonly_summary_rounded <-
  outerloop_model_tuned_var_imp_jk_test_AUC_withonly_summary %>%
  mutate(across(where(is.character), ~ as.factor(.x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x >= 0.001, round(.x, 3), .x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x < 0.001 & .x != 0, "<0.001", .x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x == 0, "0", .x)))

# save output as a csv file
write.csv(
  outerloop_model_tuned_var_imp_jk_test_AUC_withonly_summary_rounded,
  file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/outerloop_model_tuned_var_imp_jk_test_AUC_withonly_summary_rounded.csv",
  row.names = FALSE
)

# show same data frame, but this time with just the most relevant/main statistics (mean, standard deviation, and min/max) from testing each variable in isolation, then save output as a csv file
outerloop_model_tuned_var_imp_jk_test_AUC_withonly_main_stats <-
  outerloop_model_tuned_var_imp_jk_test_AUC_withonly_summary %>% dplyr::select(
    Variable,
    Variable_description,
    mean_Test_AUC_withonly,
    sd_Test_AUC_withonly,
    min_Test_AUC_withonly,
    max_Test_AUC_withonly
  )
outerloop_model_tuned_var_imp_jk_test_AUC_withonly_main_stats

write.csv(
  outerloop_model_tuned_var_imp_jk_test_AUC_withonly_main_stats,
  file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/outerloop_model_tuned_var_imp_jk_test_AUC_withonly_main_stats.csv",
  row.names = FALSE
)

# make version of table with values rounded to 3 digits
outerloop_model_tuned_var_imp_jk_test_AUC_withonly_main_stats_rounded <-
  outerloop_model_tuned_var_imp_jk_test_AUC_withonly_main_stats %>%
  mutate(across(where(is.character), ~ as.factor(.x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x >= 0.001, round(.x, 3), .x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x < 0.001 & .x != 0, "<0.001", .x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x == 0, "0", .x)))

# save output as a csv file
write.csv(
  outerloop_model_tuned_var_imp_jk_test_AUC_withonly_main_stats_rounded,
  file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/outerloop_model_tuned_var_imp_jk_test_AUC_withonly_main_stats_rounded.csv",
  row.names = FALSE
)

# adjust/transform summarized and individual data to prepare for plotting
outerloop_model_tuned_var_imp_jk_test_AUC_withonly_main_stats_plot_data <-
  outerloop_model_tuned_var_imp_jk_test_AUC_withonly_main_stats %>%
  arrange(mean_Test_AUC_withonly)
outerloop_model_tuned_var_imp_jk_test_AUC_withonly_joined_df_plot_data <-
  outerloop_model_tuned_var_imp_jk_test_AUC_withonly_joined_df %>%
  arrange(Test_AUC_withonly)

outerloop_model_tuned_var_imp_jk_test_AUC_withonly_main_stats_plot_data$Variable <-
  factor(
    outerloop_model_tuned_var_imp_jk_test_AUC_withonly_main_stats_plot_data$Variable,
    levels = outerloop_model_tuned_var_imp_jk_test_AUC_withonly_main_stats_plot_data$Variable
  )
outerloop_model_tuned_var_imp_jk_test_AUC_withonly_joined_df_plot_data$Variable <-
  factor(
    outerloop_model_tuned_var_imp_jk_test_AUC_withonly_joined_df_plot_data$Variable,
    levels = outerloop_model_tuned_var_imp_jk_test_AUC_withonly_joined_df_plot_data$Variable %>% unique()
  )

outerloop_model_tuned_var_imp_jk_test_AUC_withonly_main_stats_plot_data$Variable_description <-
  factor(
    outerloop_model_tuned_var_imp_jk_test_AUC_withonly_main_stats_plot_data$Variable_description,
    levels = outerloop_model_tuned_var_imp_jk_test_AUC_withonly_main_stats_plot_data$Variable_description
  )
outerloop_model_tuned_var_imp_jk_test_AUC_withonly_joined_df_plot_data$Variable_description <-
  factor(
    outerloop_model_tuned_var_imp_jk_test_AUC_withonly_joined_df_plot_data$Variable_description,
    levels = outerloop_model_tuned_var_imp_jk_test_AUC_withonly_joined_df_plot_data$Variable_description %>% unique()
  )

outerloop_model_tuned_var_imp_jk_test_AUC_withonly_main_stats_plot_data <-
  outerloop_model_tuned_var_imp_jk_test_AUC_withonly_main_stats_plot_data %>%
  group_by(Variable_description)
outerloop_model_tuned_var_imp_jk_test_AUC_withonly_joined_df_plot_data <-
  outerloop_model_tuned_var_imp_jk_test_AUC_withonly_joined_df_plot_data %>%
  group_by(Variable_description)

# plot bar graph showing mean univariate jackknife importance with individual sample data points added for each variable, then save the image output
(outerloop_model_tuned_var_imp_jk_test_AUC_withonly_bargraph <-
    {
      invisible(
        tmp_plot <-
          ggplot(
            data = outerloop_model_tuned_var_imp_jk_test_AUC_withonly_main_stats_plot_data,
            aes(x = Variable_description, y = mean_Test_AUC_withonly)
          ) +
          labs(x = "Predictors",
               y = "Mean jackknife importance (test AUC)") +
          coord_flip() +
          geom_bar(
            position = "dodge",
            stat = "identity",
            fill = alpha(custom_color_pal_1$green2_zomp, 0.4),
            color = alpha(custom_color_pal_1$gray_paynes, 0.9),
            linewidth = 0.6,
            width = 0.8
          ) +
          geom_point(
            data = (
              outerloop_model_tuned_var_imp_jk_test_AUC_withonly_joined_df_plot_data %>% filter(max(Test_AUC_withonly) > 0)
            ),
            aes(x = Variable_description, y = Test_AUC_withonly),
            size = 2.3,
            fill = alpha(custom_color_pal_1$yellow_saffron, 0.8),
            color = alpha(custom_color_pal_1$orange1_crayola, 0.9),
            shape = 23,
            stroke = 1.1
          ) +  # add individual data points (n = 5) for variables with at least one test AUC sample value greater than zero
          geom_hline(
            yintercept = (
              ref_line <- (
                outerloop_cv_output_stats_summary %>%
                  filter(performance_metric == "AUC") %>%
                  dplyr::select(mean)
              )$mean
            ),
            color = custom_color_pal_1$red_imperial,
            size = 1.1,
            alpha = 0.8
          ) +  # line showing reference mean test AUC (from estimated generalization error results)
          geom_text(
            aes(
              0,
              (
                outerloop_cv_output_stats_summary %>%
                  filter(performance_metric == "AUC") %>%
                  dplyr::select(mean)
              )$mean,
              label = "Reference test AUC",
              vjust = 1.5,
              hjust = -0.3,
              angle = 90
            ),
            size = 7,
            color = custom_color_pal_1$red_imperial
          ) +
          theme_bw() +
          theme(axis.title.x = element_text(
            color = "#444444",
            size = 24,
            margin = margin(20, 20, 20, 20)
          )) +
          theme(axis.title.y = element_text(
            color = "#444444",
            size = 24,
            margin = margin(20, 20, 20, 20)
          )) +
          theme(axis.text = element_text(
            color = "#444444", size = 16
          ))
      )
      tmp_plot + scale_y_continuous(
        breaks = sort(c(
          ggplot_build(tmp_plot)$layout$panel_params[[1]]$x$breaks,
          ref_line
        )),
        labels = ~ ifelse(.x == ref_line, round(.x, 3), prettyNum(.x))
      )
    })
rm(tmp_plot)

ggsave(
  file = "outerloop_model_tuned_var_imp_jk_test_AUC_withonly_bargraph.pdf",
  outerloop_model_tuned_var_imp_jk_test_AUC_withonly_bargraph,
  device = "pdf",
  path = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil",
  width = 20,
  height = 28,
  units = "in"
)
ggsave(
  file = "outerloop_model_tuned_var_imp_jk_test_AUC_withonly_bargraph.png",
  outerloop_model_tuned_var_imp_jk_test_AUC_withonly_bargraph,
  device = "png",
  path = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil",
  width = 20,
  height = 28,
  units = "in",
  dpi = 400
)

# save all objects in R environment
save.image(
  "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil.RData"
)


## Correlation analysis ----

# generate a SWD object with background sites randomly sampled from 1991-2020 temporally-averaged set of environmental raster layers, then save the environmental background SWD object
bio_env_c_bg_sample_size <- 12000

set.seed(25)
bio_env_c_bg <- randomPoints(bio_env_c,
                             n = bio_env_c_bg_sample_size)

bio_env_c_bg <- prepareSWD(species = "bg_sites",
                           a = bio_env_c_bg,
                           env = bio_env_c)

swd2csv(bio_env_c_bg, file_name = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/bio_env_c_bg.csv")

bio_env_c_bg_full_var_names <- bio_env_c_bg
names(bio_env_c_bg_full_var_names@data) <-
  bio_env_layer_names$variable_name

swd2csv(bio_env_c_bg_full_var_names, file_name = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/bio_env_c_bg_full_var_names.csv")

# plot heat map showing degree of correlation between environmental variables; save the image output
bio_env_c_bg_corvar_plot <- plotCor(bio_env_c_bg_full_var_names,
                                    method = "pearson",
                                    cor_th = 0.9)

ggsave(
  file = "bio_env_c_bg_corvar_plot.pdf",
  bio_env_c_bg_corvar_plot,
  device = "pdf",
  path = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil",
  width = 35,
  height = 35,
  units = "in"
)
ggsave(
  file = "bio_env_c_bg_corvar_plot.png",
  bio_env_c_bg_corvar_plot,
  device = "png",
  path = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil",
  width = 35,
  height = 35,
  units = "in",
  dpi = 360
)

# get data from correlation heat map; then save output as a csv file
bio_env_c_bg_corvar_plot_data <-
  bio_env_c_bg_corvar_plot$data

write.csv(bio_env_c_bg_corvar_plot_data,
          file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/bio_env_c_bg_corvar_plot_data.csv",
          row.names = FALSE)

# print pairs of correlated variables according to given method and correlation threshold; save output as a csv file
bio_env_c_bg_corvar <- corVar(
  bio_env_c_bg_full_var_names,
  method = "pearson",
  cor_th = 0.7,
  order = TRUE,
  remove_diagonal = TRUE
)
bio_env_c_bg_corvar

write.csv(bio_env_c_bg_corvar,
          file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/bio_env_c_bg_corvar.csv",
          row.names = FALSE)

# get number of pairs of significantly correlated variables (above given threshold); save output as a csv file
bio_env_c_bg_num_cor_var_pairs <- nrow(bio_env_c_bg_corvar)
print(str_c(
  "Number of correlated variable pairs: ",
  bio_env_c_bg_num_cor_var_pairs
))

write.csv(bio_env_c_bg_num_cor_var_pairs,
          file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/bio_env_c_bg_num_cor_var_pairs.csv",
          row.names = FALSE)


## VIF testing for multicollinearity ----

# calculate VIF values for all variables; environmental values sampled from a 1991-2020 temporally-averaged environmental raster set

bio_env_c_bg_vif <-
  vif(bio_env_c_bg@data,
      (if (p_version(usdm) >= "2.1-6") {
        size = nrow(bio_env_c_bg@data)
      } else{
        maxobservations = nrow(bio_env_c_bg@data)
      })) %>%
  tibble() %>%
  arrange(VIF) %>%
  left_join(
    dplyr::select(bio_env_layer_names, variable_name, variable_abbreviation),
    by = join_by("Variables" == "variable_abbreviation")
  ) %>%
  relocate(Variable_description = variable_name, .after = Variables) %>% rename(Variable = Variables)
bio_env_c_bg_vif

# save output as a csv file
write.csv(bio_env_c_bg_vif,
          file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/bio_env_c_bg_vif.csv",
          row.names = FALSE)

# make version of table with values rounded to 3 digits
bio_env_c_bg_vif_rounded <-
  bio_env_c_bg_vif %>%
  mutate(across(where(is.character), ~ as.factor(.x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x >= 0.001, round(.x, 3), .x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x < 0.001 & .x != 0, "<0.001", .x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x == 0, "0", .x)))

# calculate VIF for all variables, exclude one with highest VIF (greater than given threshold), then repeat procedure until no variables with VIF greater than given threshold remain
bio_env_c_bg_vifstep <-
  vifstep(bio_env_c_bg@data,
          (if (p_version(usdm) >= "2.1-6") {
            size = nrow(bio_env_c_bg@data)
          } else{
            maxobservations = nrow(bio_env_c_bg@data)
          }),
          th = 10)
bio_env_c_bg_vifstep
bio_env_c_bg_vifstep_excluded_vars <-
  bio_env_c_bg_vifstep@excluded
bio_env_c_bg_vifstep_corvar <-
  bio_env_c_bg_vifstep@corMatrix
bio_env_c_bg_vifstep_vif <-
  bio_env_c_bg_vifstep@results %>% arrange(VIF)

# save outputs
sink(
  "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/bio_env_c_bg_vifstep.txt"
)
print(bio_env_c_bg_vifstep)
sink()

write.csv(bio_env_c_bg_vifstep_excluded_vars,
          file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/bio_env_c_bg_vifstep_excluded_vars.csv",
          row.names = FALSE)

write.csv(bio_env_c_bg_vifstep_corvar,
          file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/bio_env_c_bg_vifstep_corvar.csv",
          row.names = FALSE)

write.csv(bio_env_c_bg_vifstep_vif %>% arrange(VIF),
          file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/bio_env_c_bg_vifstep_vif.csv",
          row.names = FALSE)

# find pair of variables that has significant linear correlation (greater than given threshold). exclude one from pair with greater VIF, then repeat procedure until no variables with correlation coefficient greater than given threshold remain
bio_env_c_bg_vifcor <-
  vifcor(bio_env_c_bg@data,
         (if (p_version(usdm) >= "2.1-6") {
           size = nrow(bio_env_c_bg@data)
         } else{
           maxobservations = nrow(bio_env_c_bg@data)
         }),
         th = 0.7)
bio_env_c_bg_vifcor
bio_env_c_bg_vifcor_excluded_vars <-
  bio_env_c_bg_vifstep@excluded
bio_env_c_bg_vifcor_corvar <-
  bio_env_c_bg_vifstep@corMatrix
bio_env_c_bg_vifcor_vif <-
  bio_env_c_bg_vifstep@results %>% arrange(VIF)

# save outputs
sink(
  "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/bio_env_c_bg_vifcor.txt"
)
print(bio_env_c_bg_vifcor)
sink()

write.csv(bio_env_c_bg_vifcor_excluded_vars,
          file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/bio_env_c_bg_vifcor_excluded_vars.csv",
          row.names = FALSE)

write.csv(bio_env_c_bg_vifcor_corvar,
          file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/bio_env_c_bg_vifcor_corvar.csv",
          row.names = FALSE)

write.csv(bio_env_c_bg_vifcor_vif %>% arrange(VIF),
          file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/bio_env_c_bg_vifcor_vif.csv",
          row.names = FALSE)


### VIF with permutation importance ----

# combine VIF and Maxent variable permutation importance results in data frame, sorted by VIF
outerloop_model_tuned_vif_and_var_imp_sorted_by_vif <-
  left_join(
    bio_env_c_bg_vif,
    outerloop_model_tuned_var_imp_main_stats,
    by = join_by(
      "Variable" == "Variable",
      "Variable_description" == "Variable_description"
    )
  ) %>%
  arrange(VIF, desc(mean_Permutation_importance))
outerloop_model_tuned_vif_and_var_imp_sorted_by_vif

# save output as a csv file
write.csv(
  outerloop_model_tuned_vif_and_var_imp_sorted_by_vif,
  file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/outerloop_model_tuned_vif_and_var_imp_sorted_by_vif.csv",
  row.names = FALSE
)

# make version of table with values rounded to 3 digits
outerloop_model_tuned_vif_and_var_imp_sorted_by_vif_rounded <-
  outerloop_model_tuned_vif_and_var_imp_sorted_by_vif %>%
  mutate(across(where(is.character), ~ as.factor(.x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x >= 0.001, round(.x, 3), .x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x < 0.001 & .x != 0, "<0.001", .x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x == 0, "0", .x)))

# save output as a csv file
write.csv(
  outerloop_model_tuned_vif_and_var_imp_sorted_by_vif_rounded,
  file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/outerloop_model_tuned_vif_and_var_imp_sorted_by_vif_rounded.csv",
  row.names = FALSE
)

# combine VIF and Maxent variable permutation importance results in data frame, sorted by mean permutation importance
outerloop_model_tuned_vif_and_var_imp_sorted_by_var_imp <-
  left_join(
    bio_env_c_bg_vif,
    outerloop_model_tuned_var_imp_main_stats,
    by = join_by(
      "Variable" == "Variable",
      "Variable_description" == "Variable_description"
    )
  ) %>%
  arrange(desc(mean_Permutation_importance), VIF)
outerloop_model_tuned_vif_and_var_imp_sorted_by_var_imp

# save output as a csv file
write.csv(
  outerloop_model_tuned_vif_and_var_imp_sorted_by_var_imp,
  file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/outerloop_model_tuned_vif_and_var_imp_sorted_by_var_imp.csv",
  row.names = FALSE
)

# make version of table with values rounded to 3 digits
outerloop_model_tuned_vif_and_var_imp_sorted_by_var_imp_rounded <-
  outerloop_model_tuned_vif_and_var_imp_sorted_by_var_imp %>%
  mutate(across(where(is.character), ~ as.factor(.x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x >= 0.001, round(.x, 3), .x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x < 0.001 & .x != 0, "<0.001", .x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x == 0, "0", .x)))

# save output as a csv file
write.csv(
  outerloop_model_tuned_vif_and_var_imp_sorted_by_var_imp_rounded,
  file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/outerloop_model_tuned_vif_and_var_imp_sorted_by_var_imp_rounded.csv",
  row.names = FALSE
)


### VIF with jackknife variable importance ----

# combine VIF and univariate jackknife variable importance results in data frame, sorted by VIF
outerloop_model_tuned_vif_and_var_imp_jk_sorted_by_vif <-
  left_join(
    bio_env_c_bg_vif,
    outerloop_model_tuned_var_imp_jk_test_AUC_withonly_main_stats,
    by = join_by(
      "Variable" == "Variable",
      "Variable_description" == "Variable_description"
    )
  ) %>%
  arrange(VIF, desc(mean_Test_AUC_withonly))
outerloop_model_tuned_vif_and_var_imp_jk_sorted_by_vif

# save output as a csv file
write.csv(
  outerloop_model_tuned_vif_and_var_imp_jk_sorted_by_vif,
  file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/outerloop_model_tuned_vif_and_var_imp_jk_sorted_by_vif.csv",
  row.names = FALSE
)

# make version of table with values rounded to 3 digits
outerloop_model_tuned_vif_and_var_imp_jk_sorted_by_vif_rounded <-
  outerloop_model_tuned_vif_and_var_imp_jk_sorted_by_vif %>%
  mutate(across(where(is.character), ~ as.factor(.x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x >= 0.001, round(.x, 3), .x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x < 0.001 & .x != 0, "<0.001", .x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x == 0, "0", .x)))

# save output as a csv file
write.csv(
  outerloop_model_tuned_vif_and_var_imp_jk_sorted_by_vif_rounded,
  file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/outerloop_model_tuned_vif_and_var_imp_jk_sorted_by_vif_rounded.csv",
  row.names = FALSE
)

# combine VIF and univariate jackknife variable importance results in data frame, sorted by mean univariate jackknife variable importance test AUC
outerloop_model_tuned_vif_and_var_imp_jk_sorted_by_var_imp_jk <-
  left_join(
    bio_env_c_bg_vif,
    outerloop_model_tuned_var_imp_jk_test_AUC_withonly_main_stats,
    by = join_by(
      "Variable" == "Variable",
      "Variable_description" == "Variable_description"
    )
  ) %>%
  arrange(desc(mean_Test_AUC_withonly), VIF)
outerloop_model_tuned_vif_and_var_imp_jk_sorted_by_var_imp_jk

# save output as a csv file
write.csv(
  outerloop_model_tuned_vif_and_var_imp_jk_sorted_by_var_imp_jk,
  file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/outerloop_model_tuned_vif_and_var_imp_sorted_by_var_imp.csv",
  row.names = FALSE
)

# make version of table with values rounded to 3 digits
outerloop_model_tuned_vif_and_var_imp_jk_sorted_by_var_imp_jk_rounded <-
  outerloop_model_tuned_vif_and_var_imp_jk_sorted_by_var_imp_jk %>%
  mutate(across(where(is.character), ~ as.factor(.x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x >= 0.001, round(.x, 3), .x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x < 0.001 & .x != 0, "<0.001", .x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x == 0, "0", .x)))

# save output as a csv file
write.csv(
  outerloop_model_tuned_vif_and_var_imp_jk_sorted_by_var_imp_jk_rounded,
  file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/outerloop_model_tuned_vif_and_var_imp_jk_sorted_by_var_imp_jk_rounded.csv",
  row.names = FALSE
)

### VIF with both variable importance measures ----

# combine VIF results with those from both variable importance measures (Maxent variable permutation importance and univariate jackknife variable importance) in data frame, sorted by VIF
outerloop_model_tuned_vif_and_both_var_imp_types_sorted_by_vif <-
  left_join(
    outerloop_model_tuned_vif_and_var_imp_jk_sorted_by_vif,
    outerloop_model_tuned_vif_and_var_imp_sorted_by_vif,
    by = join_by(
      "Variable" == "Variable",
      "Variable_description" == "Variable_description",
      "VIF" == "VIF"
    )
  ) %>%
  arrange(VIF,
          desc(mean_Test_AUC_withonly),
          desc(mean_Permutation_importance))
outerloop_model_tuned_vif_and_both_var_imp_types_sorted_by_vif

# save output as a csv file
write.csv(
  outerloop_model_tuned_vif_and_both_var_imp_types_sorted_by_vif,
  file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/outerloop_model_tuned_vif_and_both_var_imp_types_sorted_by_vif.csv",
  row.names = FALSE
)

# make version of table with values rounded to 3 digits
outerloop_model_tuned_vif_and_both_var_imp_types_sorted_by_vif_rounded <-
  outerloop_model_tuned_vif_and_both_var_imp_types_sorted_by_vif %>%
  mutate(across(where(is.character), ~ as.factor(.x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x >= 0.001, round(.x, 3), .x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x < 0.001 & .x != 0, "<0.001", .x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x == 0, "0", .x)))

# save output as a csv file
write.csv(
  outerloop_model_tuned_vif_and_both_var_imp_types_sorted_by_vif_rounded,
  file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/outerloop_model_tuned_vif_and_both_var_imp_types_sorted_by_vif_rounded.csv",
  row.names = FALSE
)

# make version of table combining VIF with only the mean values of both variable importance measures, sorted by VIF
outerloop_model_tuned_vif_and_both_var_imp_types_sorted_by_vif_means <-
  outerloop_model_tuned_vif_and_both_var_imp_types_sorted_by_vif %>%
  dplyr::select(
    Variable,
    Variable_description,
    VIF,
    mean_Test_AUC_withonly,
    mean_Permutation_importance
  )
outerloop_model_tuned_vif_and_both_var_imp_types_sorted_by_vif_means

# save output as a csv file
write.csv(
  outerloop_model_tuned_vif_and_both_var_imp_types_sorted_by_vif_means,
  file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/outerloop_model_tuned_vif_and_both_var_imp_types_sorted_by_vif_means.csv",
  row.names = FALSE
)

# make version of table with values rounded to 3 digits
outerloop_model_tuned_vif_and_both_var_imp_types_sorted_by_vif_means_rounded <-
  outerloop_model_tuned_vif_and_both_var_imp_types_sorted_by_vif_means %>%
  mutate(across(where(is.character), ~ as.factor(.x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x >= 0.001, round(.x, 3), .x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x < 0.001 & .x != 0, "<0.001", .x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x == 0, "0", .x)))

# save output as a csv file
write.csv(
  outerloop_model_tuned_vif_and_both_var_imp_types_sorted_by_vif_means_rounded,
  file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/outerloop_model_tuned_vif_and_both_var_imp_types_sorted_by_vif_means_rounded.csv",
  row.names = FALSE
)

# combine VIF results with those from both variable importance measures (Maxent variable permutation importance and univariate jackknife variable importance) in data frame, sorted by mean univariate jackknife variable importance test AUC
outerloop_model_tuned_vif_and_both_var_imp_types_sorted_by_var_imp_jk <-
  left_join(
    outerloop_model_tuned_vif_and_var_imp_jk_sorted_by_vif,
    outerloop_model_tuned_vif_and_var_imp_sorted_by_vif,
    by = join_by(
      "Variable" == "Variable",
      "Variable_description" == "Variable_description",
      "VIF" == "VIF"
    )
  ) %>%
  arrange(desc(mean_Test_AUC_withonly),
          desc(mean_Permutation_importance),
          VIF)
outerloop_model_tuned_vif_and_both_var_imp_types_sorted_by_var_imp_jk

# save output as a csv file
write.csv(
  outerloop_model_tuned_vif_and_both_var_imp_types_sorted_by_var_imp_jk,
  file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/outerloop_model_tuned_vif_and_both_var_imp_types_sorted_by_var_imp_jk.csv",
  row.names = FALSE
)

# make version of table with values rounded to 3 digits
outerloop_model_tuned_vif_and_both_var_imp_types_sorted_by_var_imp_jk_rounded <-
  outerloop_model_tuned_vif_and_both_var_imp_types_sorted_by_var_imp_jk %>%
  mutate(across(where(is.character), ~ as.factor(.x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x >= 0.001, round(.x, 3), .x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x < 0.001 & .x != 0, "<0.001", .x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x == 0, "0", .x)))

# save output as a csv file
write.csv(
  outerloop_model_tuned_vif_and_both_var_imp_types_sorted_by_var_imp_jk_rounded,
  file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/outerloop_model_tuned_vif_and_both_var_imp_types_sorted_by_var_imp_jk_rounded.csv",
  row.names = FALSE
)

# make version of table combining VIF with only the mean values of both variable importance measures, sorted by mean univariate jackknife variable importance test AUC
outerloop_model_tuned_vif_and_both_var_imp_types_sorted_by_var_imp_jk_means <-
  outerloop_model_tuned_vif_and_both_var_imp_types_sorted_by_var_imp_jk %>%
  dplyr::select(
    Variable,
    Variable_description,
    VIF,
    mean_Test_AUC_withonly,
    mean_Permutation_importance
  )
outerloop_model_tuned_vif_and_both_var_imp_types_sorted_by_var_imp_jk_means

# save output as a csv file
write.csv(
  outerloop_model_tuned_vif_and_both_var_imp_types_sorted_by_var_imp_jk_means,
  file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/outerloop_model_tuned_vif_and_both_var_imp_types_sorted_by_var_imp_jk_means.csv",
  row.names = FALSE
)

# make version of table with values rounded to 3 digits
outerloop_model_tuned_vif_and_both_var_imp_types_sorted_by_var_imp_jk_means_rounded <-
  outerloop_model_tuned_vif_and_both_var_imp_types_sorted_by_var_imp_jk_means %>%
  mutate(across(where(is.character), ~ as.factor(.x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x >= 0.001, round(.x, 3), .x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x < 0.001 & .x != 0, "<0.001", .x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x == 0, "0", .x)))

# save output as a csv file
write.csv(
  outerloop_model_tuned_vif_and_both_var_imp_types_sorted_by_var_imp_jk_means_rounded,
  file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/outerloop_model_tuned_vif_and_both_var_imp_types_sorted_by_var_imp_jk_means_rounded.csv",
  row.names = FALSE
)

# combine VIF results with those from both variable importance measures (Maxent variable permutation importance and univariate jackknife variable importance) in data frame, sorted by mean permutation importance
outerloop_model_tuned_vif_and_both_var_imp_types_sorted_by_var_imp <-
  left_join(
    outerloop_model_tuned_vif_and_var_imp_sorted_by_vif,
    outerloop_model_tuned_vif_and_var_imp_jk_sorted_by_vif,
    by = join_by(
      "Variable" == "Variable",
      "Variable_description" == "Variable_description",
      "VIF" == "VIF"
    )
  ) %>%
  arrange(desc(mean_Permutation_importance),
          desc(mean_Test_AUC_withonly),
          VIF)
outerloop_model_tuned_vif_and_both_var_imp_types_sorted_by_var_imp

# save output as a csv file
write.csv(
  outerloop_model_tuned_vif_and_both_var_imp_types_sorted_by_var_imp,
  file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/outerloop_model_tuned_vif_and_both_var_imp_types_sorted_by_var_imp.csv",
  row.names = FALSE
)

# make version of table with values rounded to 3 digits
outerloop_model_tuned_vif_and_both_var_imp_types_sorted_by_var_imp_rounded <-
  outerloop_model_tuned_vif_and_both_var_imp_types_sorted_by_var_imp %>%
  mutate(across(where(is.character), ~ as.factor(.x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x >= 0.001, round(.x, 3), .x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x < 0.001 & .x != 0, "<0.001", .x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x == 0, "0", .x)))

# save output as a csv file
write.csv(
  outerloop_model_tuned_vif_and_both_var_imp_types_sorted_by_var_imp_rounded,
  file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/outerloop_model_tuned_vif_and_both_var_imp_types_sorted_by_var_imp_rounded.csv",
  row.names = FALSE
)

# make version of table combining VIF with only the mean values of both variable importance measures, sorted by mean permutation importance
outerloop_model_tuned_vif_and_both_var_imp_types_sorted_by_var_imp_means <-
  outerloop_model_tuned_vif_and_both_var_imp_types_sorted_by_var_imp %>%
  dplyr::select(
    Variable,
    Variable_description,
    VIF,
    mean_Permutation_importance,
    mean_Test_AUC_withonly
  )
outerloop_model_tuned_vif_and_both_var_imp_types_sorted_by_var_imp_means

# save output as a csv file
write.csv(
  outerloop_model_tuned_vif_and_both_var_imp_types_sorted_by_var_imp_means,
  file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/outerloop_model_tuned_vif_and_both_var_imp_types_sorted_by_var_imp_means.csv",
  row.names = FALSE
)

# make version of table with values rounded to 3 digits
outerloop_model_tuned_vif_and_both_var_imp_types_sorted_by_var_imp_means_rounded <-
  outerloop_model_tuned_vif_and_both_var_imp_types_sorted_by_var_imp_means %>%
  mutate(across(where(is.character), ~ as.factor(.x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x >= 0.001, round(.x, 3), .x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x < 0.001 & .x != 0, "<0.001", .x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x == 0, "0", .x)))

# save output as a csv file
write.csv(
  outerloop_model_tuned_vif_and_both_var_imp_types_sorted_by_var_imp_means_rounded,
  file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/outerloop_model_tuned_vif_and_both_var_imp_types_sorted_by_var_imp_means_rounded.csv",
  row.names = FALSE
)


### VIFstep subset with both variable importance measures ----

# combine stepwise VIF results with those from both variable importance measures (Maxent variable permutation importance and univariate jackknife variable importance) in data frame, sorted by VIF; exclude multicollinear variables (VIF > 10) per results from vifstep algorithm (which involved stepwise VIF calculations per variable)
outerloop_model_tuned_vifstep_subset_and_both_var_imp_types_sorted_by_vif <-
  left_join(
    outerloop_model_tuned_vif_and_both_var_imp_types_sorted_by_vif %>% dplyr::select(-VIF) %>% filter(!(
      Variable %in% bio_env_c_bg_vifstep_excluded_vars
    )),
    bio_env_c_bg_vifstep_vif,
    join_by("Variable" == "Variables")
  ) %>% relocate(VIF, .after = Variable_description) %>%
  arrange(VIF,
          desc(mean_Test_AUC_withonly),
          desc(mean_Permutation_importance))
outerloop_model_tuned_vifstep_subset_and_both_var_imp_types_sorted_by_vif

# save output as a csv file
write.csv(
  outerloop_model_tuned_vifstep_subset_and_both_var_imp_types_sorted_by_vif,
  file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/outerloop_model_tuned_vifstep_subset_and_both_var_imp_types_sorted_by_vif.csv",
  row.names = FALSE
)

# make version of table with values rounded to 3 digits
outerloop_model_tuned_vifstep_subset_and_both_var_imp_types_sorted_by_vif_rounded <-
  outerloop_model_tuned_vifstep_subset_and_both_var_imp_types_sorted_by_vif %>%
  mutate(across(where(is.character), ~ as.factor(.x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x >= 0.001, round(.x, 3), .x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x < 0.001 & .x != 0, "<0.001", .x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x == 0, "0", .x)))

# save output as a csv file
write.csv(
  outerloop_model_tuned_vifstep_subset_and_both_var_imp_types_sorted_by_vif_rounded,
  file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/outerloop_model_tuned_vifstep_subset_and_both_var_imp_types_sorted_by_vif_rounded.csv",
  row.names = FALSE
)

# make version of table combining stepwise VIF with only the mean values of both variable importance measures, sorted by VIF
outerloop_model_tuned_vifstep_subset_and_both_var_imp_types_sorted_by_vif_means <-
  outerloop_model_tuned_vifstep_subset_and_both_var_imp_types_sorted_by_vif %>%
  dplyr::select(
    Variable,
    Variable_description,
    VIF,
    mean_Test_AUC_withonly,
    mean_Permutation_importance
  )
outerloop_model_tuned_vifstep_subset_and_both_var_imp_types_sorted_by_vif_means

# save output as a csv file
write.csv(
  outerloop_model_tuned_vifstep_subset_and_both_var_imp_types_sorted_by_vif_means,
  file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/outerloop_model_tuned_vifstep_subset_and_both_var_imp_types_sorted_by_vif_means.csv",
  row.names = FALSE
)

# make version of table with values rounded to 3 digits
outerloop_model_tuned_vifstep_subset_and_both_var_imp_types_sorted_by_vif_means_rounded <-
  outerloop_model_tuned_vifstep_subset_and_both_var_imp_types_sorted_by_vif_means %>%
  mutate(across(where(is.character), ~ as.factor(.x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x >= 0.001, round(.x, 3), .x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x < 0.001 & .x != 0, "<0.001", .x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x == 0, "0", .x)))

# save output as a csv file
write.csv(
  outerloop_model_tuned_vifstep_subset_and_both_var_imp_types_sorted_by_vif_means_rounded,
  file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/outerloop_model_tuned_vifstep_subset_and_both_var_imp_types_sorted_by_vif_means_rounded.csv",
  row.names = FALSE
)

# combine stepwise VIF results with those from both variable importance measures (Maxent variable permutation importance and univariate jackknife variable importance) in data frame, sorted by mean univariate jackknife variable importance test AUC; exclude multicollinear variables (VIF > 10) per results from vifstep algorithm (which involved stepwise VIF calculations per variable)
outerloop_model_tuned_vifstep_subset_and_both_var_imp_types_sorted_by_var_imp_jk <-
  left_join(
    outerloop_model_tuned_vif_and_both_var_imp_types_sorted_by_var_imp_jk %>% dplyr::select(-VIF) %>% filter(!(
      Variable %in% bio_env_c_bg_vifstep_excluded_vars
    )),
    bio_env_c_bg_vifstep_vif,
    join_by("Variable" == "Variables")
  ) %>% relocate(VIF, .after = Variable_description) %>%
  arrange(desc(mean_Test_AUC_withonly),
          desc(mean_Permutation_importance),
          VIF)
outerloop_model_tuned_vifstep_subset_and_both_var_imp_types_sorted_by_var_imp_jk

# save output as a csv file
write.csv(
  outerloop_model_tuned_vifstep_subset_and_both_var_imp_types_sorted_by_var_imp_jk,
  file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/outerloop_model_tuned_vifstep_subset_and_both_var_imp_types_sorted_by_var_imp_jk.csv",
  row.names = FALSE
)

# make version of table with values rounded to 3 digits
outerloop_model_tuned_vifstep_subset_and_both_var_imp_types_sorted_by_var_imp_jk_rounded <-
  outerloop_model_tuned_vifstep_subset_and_both_var_imp_types_sorted_by_var_imp_jk %>%
  mutate(across(where(is.character), ~ as.factor(.x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x >= 0.001, round(.x, 3), .x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x < 0.001 & .x != 0, "<0.001", .x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x == 0, "0", .x)))

# save output as a csv file
write.csv(
  outerloop_model_tuned_vifstep_subset_and_both_var_imp_types_sorted_by_var_imp_jk_rounded,
  file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/outerloop_model_tuned_vifstep_subset_and_both_var_imp_types_sorted_by_var_imp_jk_rounded.csv",
  row.names = FALSE
)

# make version of table combining stepwise VIF with only the mean values of both variable importance measures, sorted by mean univariate jackknife variable importance test AUC
outerloop_model_tuned_vifstep_subset_and_both_var_imp_types_sorted_by_var_imp_jk_means <-
  outerloop_model_tuned_vifstep_subset_and_both_var_imp_types_sorted_by_var_imp_jk %>%
  dplyr::select(
    Variable,
    Variable_description,
    VIF,
    mean_Test_AUC_withonly,
    mean_Permutation_importance
  )
outerloop_model_tuned_vifstep_subset_and_both_var_imp_types_sorted_by_var_imp_jk_means

# save output as a csv file
write.csv(
  outerloop_model_tuned_vifstep_subset_and_both_var_imp_types_sorted_by_var_imp_jk_means,
  file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/outerloop_model_tuned_vifstep_subset_and_both_var_imp_types_sorted_by_var_imp_jk_means.csv",
  row.names = FALSE
)

# make version of table with values rounded to 3 digits
outerloop_model_tuned_vifstep_subset_and_both_var_imp_types_sorted_by_var_imp_jk_means_rounded <-
  outerloop_model_tuned_vifstep_subset_and_both_var_imp_types_sorted_by_var_imp_jk_means %>%
  mutate(across(where(is.character), ~ as.factor(.x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x >= 0.001, round(.x, 3), .x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x < 0.001 & .x != 0, "<0.001", .x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x == 0, "0", .x)))

# save output as a csv file
write.csv(
  outerloop_model_tuned_vifstep_subset_and_both_var_imp_types_sorted_by_var_imp_jk_means_rounded,
  file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/outerloop_model_tuned_vifstep_subset_and_both_var_imp_types_sorted_by_var_imp_jk_means_rounded.csv",
  row.names = FALSE
)

# combine stepwise VIF results with those from both variable importance measures (Maxent variable permutation importance and univariate jackknife variable importance) in data frame, sorted by mean permutation importance; exclude multicollinear variables (VIF > 10) per results from vifstep algorithm (which involved stepwise VIF calculations per variable)
outerloop_model_tuned_vifstep_subset_and_both_var_imp_types_sorted_by_var_imp <-
  left_join(
    outerloop_model_tuned_vif_and_both_var_imp_types_sorted_by_var_imp %>% dplyr::select(-VIF) %>% filter(!(
      Variable %in% bio_env_c_bg_vifstep_excluded_vars
    )),
    bio_env_c_bg_vifstep_vif,
    join_by("Variable" == "Variables")
  ) %>% relocate(VIF, .after = Variable_description) %>%
  arrange(desc(mean_Permutation_importance),
          desc(mean_Test_AUC_withonly),
          VIF)
outerloop_model_tuned_vifstep_subset_and_both_var_imp_types_sorted_by_var_imp

# save output as a csv file
write.csv(
  outerloop_model_tuned_vifstep_subset_and_both_var_imp_types_sorted_by_var_imp,
  file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/outerloop_model_tuned_vifstep_subset_and_both_var_imp_types_sorted_by_var_imp.csv",
  row.names = FALSE
)

# make version of table with values rounded to 3 digits
outerloop_model_tuned_vifstep_subset_and_both_var_imp_types_sorted_by_var_imp_rounded <-
  outerloop_model_tuned_vifstep_subset_and_both_var_imp_types_sorted_by_var_imp %>%
  mutate(across(where(is.character), ~ as.factor(.x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x >= 0.001, round(.x, 3), .x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x < 0.001 & .x != 0, "<0.001", .x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x == 0, "0", .x)))

# save output as a csv file
write.csv(
  outerloop_model_tuned_vifstep_subset_and_both_var_imp_types_sorted_by_var_imp_rounded,
  file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/outerloop_model_tuned_vifstep_subset_and_both_var_imp_types_sorted_by_var_imp_rounded.csv",
  row.names = FALSE
)

# make version of table combining stepwise VIF with only the mean values of both variable importance measures, sorted by mean permutation importance
outerloop_model_tuned_vifstep_subset_and_both_var_imp_types_sorted_by_var_imp_means <-
  outerloop_model_tuned_vifstep_subset_and_both_var_imp_types_sorted_by_var_imp %>%
  dplyr::select(
    Variable,
    Variable_description,
    VIF,
    mean_Permutation_importance,
    mean_Test_AUC_withonly
  )
outerloop_model_tuned_vifstep_subset_and_both_var_imp_types_sorted_by_var_imp_means

# save output as a csv file
write.csv(
  outerloop_model_tuned_vifstep_subset_and_both_var_imp_types_sorted_by_var_imp_means,
  file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/outerloop_model_tuned_vifstep_subset_and_both_var_imp_types_sorted_by_var_imp_means.csv",
  row.names = FALSE
)

# make version of table with values rounded to 3 digits
outerloop_model_tuned_vifstep_subset_and_both_var_imp_types_sorted_by_var_imp_means_rounded <-
  outerloop_model_tuned_vifstep_subset_and_both_var_imp_types_sorted_by_var_imp_means %>%
  mutate(across(where(is.character), ~ as.factor(.x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x >= 0.001, round(.x, 3), .x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x < 0.001 & .x != 0, "<0.001", .x))) %>%
  mutate(across(!where(is.factor),
                ~ ifelse(.x == 0, "0", .x)))

# save output as a csv file
write.csv(
  outerloop_model_tuned_vifstep_subset_and_both_var_imp_types_sorted_by_var_imp_means_rounded,
  file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/outerloop_model_tuned_vifstep_subset_and_both_var_imp_types_sorted_by_var_imp_means_rounded.csv",
  row.names = FALSE
)


# Use final model to create predictions and estimated distribution maps ----

## Generate predictions from final model and environmental raster data ----

# update list of time step values to remove those corresponding to SWD objects that had dropped all presence points due to insufficient environmental data, then remove those zero-presence SWD objects
{
  if (exists("timesteps_old")) {
    timesteps <- timesteps_old
  }
  timesteps_old <- timesteps

  if (exists("data_list_old")) {
    data_list <- data_list_old
  }
  data_list_old <- data_list
  data_list_new <- data_list

  for (i in seq_along(data_list)) {
    if (length(data_list[[i]]@pa[data_list[[i]]@pa == 1]) == 0) {
      timesteps <- timesteps[!timesteps %in% timesteps[[i]]]
      data_list_new[[i]] <- NULL
    }
  }

  data_list <- data_list_new
}

# prepare list of time steps values for temporally dynamic predictions
timesteps_all <-
  as.character(seq(1910, 2010, 10))
timesteps_no_presences <- setdiff(timesteps_all, timesteps)

# define presence localities that were included in pooled SWD object (i.e., excluding points dropped due to missing environmental data)
presence_final <- data@coords[data@pa == 1,] %>%
  rename(longitude = X,
         latitude = Y)

# define presence localities that were included in time series SWD objects (i.e., excluding points dropped due to missing environmental data)
{
  presence_final_list <- NULL

  for (i in seq_along(timesteps)) {
    presence_final_list[[i]] <-
      data_list[[i]]@coords[data_list[[i]]@pa == 1,] %>%
      rename(longitude = X,
             latitude = Y)
  }

  # rename time-specific list elements for clarity
  names(presence_final_list) <-
    paste0("presence_final_",
           timesteps)
  }

# make sf object versions for final presence coords
presence_final_sf <- st_as_sf(
  x = presence_final,
  coords = c("longitude", "latitude"),
  crs = st_crs(elevation)
)

presence_final_sf_list <- lapply(
  presence_final_list,
  st_as_sf,
  coords = c("longitude", "latitude"),
  crs = st_crs(elevation)
)


### Prediction: bio_env_c ----

#### Static ----

# generate an SDM prediction raster from the bio_env_c rasterstack object
set.seed(25)
final_prediction <-
  SDMtune::predict(
    final_model,
    data = bio_env_c,
    type = "cloglog",
    file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/final_prediction.tif",
    overwrite = TRUE
  ) %>% raster::raster()


#### Time series ----

# same as above, but this time generate a time series of SDM prediction rasters based on time-specific rasterstack objects

{
  plan(multisession, workers = num_parallel_workers)
}  # initiate parallel processing

# import raster files (environmental variables/predictors) for time steps without presence records
files_by_timestep_no_presences_list <-
  foreach(i = seq_along(timesteps_no_presences),
          .options.future = list(seed = TRUE)) %dofuture% {
            files_by_timestep_no_presences <- c(
              list.files(
                dir(
                  str_c(
                    main_path,
                    "environmental_data/NEW_environmental_data/ClimateNA/final/decadal_1911_2020"
                  ),
                  pattern = str_c("Decade_", str_sub(timesteps_no_presences[[i]], end = 3)),
                  full.names = TRUE,
                  recursive = TRUE,
                  include.dirs = TRUE
                ),
                pattern = "tif$",
                full.names = TRUE,
                recursive = TRUE
              ),
              list.files(
                c(
                  str_c(
                    main_path,
                    "environmental_data/NEW_environmental_data/EarthEnv_Consensus_Land_Cover/final"
                  ),
                  str_c(
                    main_path,
                    "environmental_data/NEW_environmental_data/HydroSHEDS_Hydrology/final"
                  ),
                  str_c(
                    main_path,
                    "environmental_data/NEW_environmental_data/SoilGrids/final_NAfilled"
                  )
                ),
                pattern = 'tif$',
                full.names = TRUE,
                recursive = TRUE
              )
            )
            files_by_timestep_no_presences <-
              setdiff(
                files_by_timestep_no_presences,
                str_subset(files_by_timestep_no_presences, "monthly")
              )
            files_by_timestep_no_presences
          }

# stack raster files (environmental variables/predictors) for time steps without presence records
bio_env_c_no_presences_list <-
  foreach(i = seq_along(timesteps_no_presences),
          .options.future = list(seed = TRUE)) %dofuture% {
            files_by_timestep_no_presences_list[[i]] %>%
              raster::stack() %>%
              raster::subset(bio_env_layer_names$variable_abbreviation) %>%
              raster::mask(as_Spatial(study_area_mask)) %>%
              raster::trim(
                filename = str_c(
                  raster_tempdir,
                  "bio_env_c_",
                  timesteps_no_presences[[i]],
                  ".grd"
                ),
                overwrite = TRUE
              )
          }

# rename time-specific list elements for clarity
names(bio_env_c_no_presences_list) <-
  paste0("bio_env_c_", timesteps_no_presences)
summary(bio_env_c_no_presences_list)

# combine raster sets from all time steps
bio_env_c_all_list <-
  c(bio_env_c_list, bio_env_c_no_presences_list)[order(names(c(
    bio_env_c_list, bio_env_c_no_presences_list
  )))]

{
  plan(sequential)
}  # end parallel processing

{
  final_prediction_list <-
    foreach(i = seq_along(timesteps_all)) %do% {
      set.seed(25)
      SDMtune::predict(
        object = final_model,
        data = bio_env_c_all_list[[i]],
        type = "cloglog",
        file = str_c(
          "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/final_prediction_",
          timesteps_all[[i]],
          ".tif"
        ),
        overwrite = TRUE
      ) %>% raster::raster()
    }

  # rename time-specific list elements for clarity
  names(final_prediction_list) <-
    paste0("final_prediction_", timesteps_all)
}


### Prediction: bio_env ----

# same as above, but this time use the full bio_env rasterstack object instead of bio_env_c, thus enabling an expanded view of the estimated distribution maps beyond the study area and allowing for more spatial context


#### Static ----

set.seed(25)
final_prediction_expanded <-
  SDMtune::predict(
    final_model,
    data = bio_env,
    type = "cloglog",
    file = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/final_prediction_expanded.tif",
    overwrite = TRUE
  ) %>% raster::raster()


#### Time series ----

# same as above, but this time generate a time series of SDM prediction rasters based on time-specific rasterstack objects

{
  plan(multisession, workers = num_parallel_workers)
}  # initiate parallel processing

# import raster files (environmental variables/predictors) for time steps without presence records
files_by_all_timesteps_list <-
  foreach(i = seq_along(timesteps_all),
          .options.future = list(seed = TRUE)) %dofuture% {
            files_by_all_timesteps <- c(
              list.files(
                dir(
                  str_c(
                    main_path,
                    "environmental_data/NEW_environmental_data/ClimateNA/final/decadal_1911_2020"
                  ),
                  pattern = str_c("Decade_", str_sub(timesteps_all[[i]], end = 3)),
                  full.names = TRUE,
                  recursive = TRUE,
                  include.dirs = TRUE
                ),
                pattern = "tif$",
                full.names = TRUE,
                recursive = TRUE
              ),
              list.files(
                c(
                  str_c(
                    main_path,
                    "environmental_data/NEW_environmental_data/EarthEnv_Consensus_Land_Cover/final"
                  ),
                  str_c(
                    main_path,
                    "environmental_data/NEW_environmental_data/HydroSHEDS_Hydrology/final"
                  ),
                  str_c(
                    main_path,
                    "environmental_data/NEW_environmental_data/SoilGrids/final_NAfilled"
                  )
                ),
                pattern = 'tif$',
                full.names = TRUE,
                recursive = TRUE
              )
            )
            files_by_all_timesteps <-
              setdiff(files_by_all_timesteps,
                      str_subset(files_by_all_timesteps, "monthly"))
            files_by_all_timesteps
          }

# stack raster files (environmental variables/predictors) for all time steps
bio_env_list <-
  foreach(i = seq_along(timesteps_all),
          .options.future = list(seed = TRUE)) %dofuture% {
            files_by_all_timesteps_list[[i]] %>%
              raster::stack() %>%
              raster::subset(bio_env_layer_names$variable_abbreviation) %>%
              raster::trim(
                filename = str_c(raster_tempdir,
                                 "bio_env_",
                                 timesteps_all[[i]],
                                 ".grd"),
                overwrite = TRUE
              )
          }

# rename time-specific list elements for clarity
names(bio_env_list) <- paste0("bio_env_", timesteps_all)
summary(bio_env_list)

{
  plan(sequential)
}  # end parallel processing

{
  final_prediction_expanded_list <-
    foreach(i = seq_along(timesteps_all)) %do% {
      set.seed(25)
      SDMtune::predict(
        object = final_model,
        data = bio_env_list[[i]],
        type = "cloglog",
        file = str_c(
          "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/final_prediction_expanded_",
          timesteps_all[[i]],
          ".tif"
        ),
        overwrite = TRUE
      ) %>% raster::raster()
    }

  # rename time-specific list elements for clarity
  names(final_prediction_expanded_list) <-
    paste0("final_prediction_expanded_",
           timesteps_all)
}

# save all objects in R environment
save.image(
  "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil.RData"
)


## Plot continuous distribution maps ----

### Continuous maps: final_prediction ----

#### Static ----

# create continuous distribution maps

study_area_plot_layer <- {
  tm_shape(study_area_mask,
           projection = st_crs(study_area_mask),
           bbox = study_area_mask) + tm_polygons(
             alpha = 0,
             border.alpha = 0.3,
             lwd = 0.7,
             border.col = "grey65"
           )
}

final_distribution_map_continuous_no_points_layer <- {
  tm_shape(final_prediction,
           raster.downsample = FALSE,
           is.master = TRUE) + tm_raster(
             # palette = viridis(
             #   256,
             #   begin = 0.3,
             #   end = 1,
             #   option = "mako",
             #   direction = -1
             # ),
             # palette = viridis(
             #   256,
             #   begin = 0.1,
             #   end = 0.9,
             #   option = "magma",
             #   direction = 1
             # ),
             palette = viridis(
               256,
               begin = 0.1,
               end = 0.9,
               option = "inferno",
               direction = 1
             ),
             # palette = tmaptools::get_brewer_pal(
             #   "Purples",
             #   n = 256,
             #   contrast = c(0.15, 0.8),
             #   plot = FALSE
             # ),
             alpha = 0.85,
             title = "Habitat\nsuitability",
             style = "cont",
             breaks = seq(0, 1, 0.1),
             legend.reverse = TRUE
           ) + tm_legend(
             position = c("left", "bottom"),
             frame = "lightgrey",
             legend.format = list(
               fun = function(x) {
                 ifelse(x %in% c(0.00, 0.25, 0.50, 0.75, 1.00), x, "")
               }
             ),
             bg.alpha = 0.8
           )
}

final_distribution_map_continuous_layer <- {
  final_distribution_map_continuous_no_points_layer + tm_shape(presence_final_sf) + tm_dots(
    size = 0.1,
    shape = 21,
    col = "white",
    alpha = 0.8,
    border.lwd = 0.6,
    border.alpha = 0.95
  )
}

(
  final_distribution_map_continuous_no_points <-
    boundaries_plot_layer + final_distribution_map_continuous_no_points_layer
)

(
  final_distribution_map_continuous <-
    boundaries_plot_layer + final_distribution_map_continuous_layer
)

# save the image output
tmap_save(
  final_distribution_map_continuous_no_points,
  filename = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/final_distribution_map_continuous_no_points.png",
  scale = 1.3
)

tmap_save(final_distribution_map_continuous,
          filename = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/final_distribution_map_continuous.png",
          scale = 1.3)


#### Time series ----

# same as above, but this time generate a time series of maps based on time-specific prediction rasters

{
  plan(multisession, workers = num_parallel_workers)
}  # initiate parallel processing

{
  final_distribution_map_continuous_no_points_layer_list <-
    foreach(i = seq_along(timesteps_all),
            .options.future = list(seed = TRUE)) %dofuture% {
              tm_shape(final_prediction_list[[i]],
                       raster.downsample = FALSE,
                       is.master = TRUE) + tm_raster(
                         palette = viridis(
                           256,
                           begin = 0.1,
                           end = 0.9,
                           option = "inferno",
                           direction = 1
                         ),
                         alpha = 0.85,
                         title = "Habitat\nsuitability",
                         style = "cont",
                         breaks = seq(0, 1, 0.1),
                         legend.reverse = TRUE
                       ) + tm_legend(
                         position = c("left", "bottom"),
                         frame = "lightgrey",
                         legend.format = list(
                           fun = function(x) {
                             ifelse(x %in% c(0.00, 0.25, 0.50, 0.75, 1.00), x, "")
                           }
                         ),
                         bg.alpha = 0.8
                       ) +
                tm_layout(panel.labels = timesteps_all[[i]])  # plot title (panel label format)
            }

  # rename time-specific list elements for clarity
  names(final_distribution_map_continuous_no_points_layer_list) <-
    paste0("final_distribution_map_continuous_no_points_layer_",
           timesteps_all)
}

{
  final_distribution_map_continuous_layer_list <-
    foreach(
      i = str_which(
        timesteps_all,
        str_c(timesteps_no_presences, collapse = "|"),
        negate = TRUE
      ),
      j = seq_along(timesteps),
      .options.future = list(seed = TRUE)
    ) %dofuture% {
      final_distribution_map_continuous_no_points_layer_list[[i]] + tm_shape(presence_final_sf_list[[j]]) + tm_dots(
        size = 0.1,
        shape = 21,
        col = "white",
        alpha = 0.8,
        border.lwd = 0.6,
        border.alpha = 0.95
      )
    }

  # rename time-specific list elements for clarity
  names(final_distribution_map_continuous_layer_list) <-
    paste0("final_distribution_map_continuous_layer_", timesteps)
}

{
  final_distribution_map_continuous_no_points_list <-
    foreach(i = seq_along(timesteps_all),
            .options.future = list(seed = TRUE)) %dofuture% {
              final_distribution_map_continuous_no_points_by_timestep <-
                {
                  boundaries_plot_layer + final_distribution_map_continuous_no_points_layer_list[[i]]
                }

              # save the image output
              tmap_save(
                final_distribution_map_continuous_no_points_by_timestep,
                filename = str_c(
                  "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/final_distribution_map_continuous_no_points_",
                  timesteps_all[[i]],
                  ".png"
                ),
                scale = 1.3
              )

              final_distribution_map_continuous_no_points_by_timestep
            }

  # rename time-specific list elements for clarity
  names(final_distribution_map_continuous_no_points_list) <-
    paste0("final_distribution_map_continuous_no_points_",
           timesteps_all)
}

{
  final_distribution_map_continuous_list <-
    foreach(i = seq_along(timesteps),
            .options.future = list(seed = TRUE)) %dofuture% {
              final_distribution_map_continuous_by_timestep <-
                {
                  boundaries_plot_layer + final_distribution_map_continuous_layer_list[[i]]
                }

              # save the image output
              tmap_save(
                final_distribution_map_continuous_by_timestep,
                filename = str_c(
                  "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/final_distribution_map_continuous_",
                  timesteps[[i]],
                  ".png"
                ),
                scale = 1.3
              )

              final_distribution_map_continuous_by_timestep
            }

  # rename time-specific list elements for clarity
  names(final_distribution_map_continuous_list) <-
    paste0("final_distribution_map_continuous_",
           timesteps)
}

{
  plan(sequential)
}  # end parallel processing


### Continuous maps: final_prediction_expanded ----

# same plotting steps as above, but this time use final_prediction_expanded (created from bio_env) instead of final_prediction (created from bio_env_c), thus enabling an expanded view of the estimated distribution maps beyond the study area and allowing for more spatial context


#### Static ----

final_distribution_expanded_map_continuous_no_points_layer <- {
  tm_shape(final_prediction_expanded,
           raster.downsample = FALSE,
           is.master = TRUE) + tm_raster(
             palette = viridis(
               256,
               begin = 0.1,
               end = 0.9,
               option = "inferno",
               direction = 1
             ),
             alpha = 0.85,
             title = "Habitat\nsuitability",
             style = "cont",
             breaks = seq(0, 1, 0.1),
             legend.reverse = TRUE
           ) + tm_legend(
             position = c("left", "bottom"),
             frame = "lightgrey",
             legend.format = list(
               fun = function(x) {
                 ifelse(x %in% c(0.00, 0.25, 0.50, 0.75, 1.00), x, "")
               }
             ),
             bg.alpha = 0.8
           )
}

final_distribution_expanded_map_continuous_layer <- {
  final_distribution_expanded_map_continuous_no_points_layer + tm_shape(presence_final_sf) + tm_dots(
    size = 0.1,
    shape = 21,
    col = "white",
    alpha = 0.8,
    border.lwd = 0.6,
    border.alpha = 0.95
  )
}

(
  final_distribution_expanded_map_continuous_no_points <-
    boundaries_plot_layer + final_distribution_expanded_map_continuous_no_points_layer
)

(
  final_distribution_expanded_map_continuous <-
    boundaries_plot_layer + final_distribution_expanded_map_continuous_layer
)

# save the image output
tmap_save(
  final_distribution_expanded_map_continuous_no_points,
  filename = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/final_distribution_expanded_map_continuous_no_points.png",
  scale = 1.3
)

tmap_save(
  final_distribution_expanded_map_continuous,
  filename = "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/final_distribution_expanded_map_continuous.png",
  scale = 1.3
)


#### Time series ----

# same as above, but this time generate a time series of maps based on time-specific prediction rasters

{
  plan(multisession, workers = num_parallel_workers)
}  # initiate parallel processing

{
  final_distribution_expanded_map_continuous_no_points_layer_list <-
    foreach(i = seq_along(timesteps_all),
            .options.future = list(seed = TRUE)) %dofuture% {
              tm_shape(
                final_prediction_expanded_list[[i]],
                raster.downsample = FALSE,
                is.master = TRUE
              ) + tm_raster(
                palette = viridis(
                  256,
                  begin = 0.1,
                  end = 0.9,
                  option = "inferno",
                  direction = 1
                ),
                alpha = 0.85,
                title = "Habitat\nsuitability",
                style = "cont",
                breaks = seq(0, 1, 0.1),
                legend.reverse = TRUE
              ) + tm_legend(
                position = c("left", "bottom"),
                frame = "lightgrey",
                legend.format = list(
                  fun = function(x) {
                    ifelse(x %in% c(0.00, 0.25, 0.50, 0.75, 1.00), x, "")
                  }
                ),
                bg.alpha = 0.8
              ) +
                tm_layout(panel.labels = timesteps_all[[i]])  # plot title (panel label format)
            }

  # rename time-specific list elements for clarity
  names(final_distribution_expanded_map_continuous_no_points_layer_list) <-
    paste0("final_distribution_map_continuous_no_points_layer_",
           timesteps_all)
}

{
  final_distribution_expanded_map_continuous_layer_list <-
    foreach(
      i = str_which(
        timesteps_all,
        str_c(timesteps_no_presences, collapse = "|"),
        negate = TRUE
      ),
      j = seq_along(timesteps),
      .options.future = list(seed = TRUE)
    ) %dofuture% {
      final_distribution_expanded_map_continuous_no_points_layer_list[[i]] + tm_shape(presence_final_sf_list[[j]]) + tm_dots(
        size = 0.1,
        shape = 21,
        col = "white",
        alpha = 0.8,
        border.lwd = 0.6,
        border.alpha = 0.95
      )
    }

  # rename time-specific list elements for clarity
  names(final_distribution_expanded_map_continuous_layer_list) <-
    paste0("final_distribution_map_continuous_layer_", timesteps)
}

{
  final_distribution_expanded_map_continuous_no_points_list <-
    foreach(i = seq_along(timesteps_all),
            .options.future = list(seed = TRUE)) %dofuture% {
              final_distribution_expanded_map_continuous_no_points_by_timestep <-
                {
                  boundaries_plot_layer + final_distribution_expanded_map_continuous_no_points_layer_list[[i]]
                }

              # save the image output
              tmap_save(
                final_distribution_expanded_map_continuous_no_points_by_timestep,
                filename = str_c(
                  "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/final_distribution_expanded_map_continuous_no_points_",
                  timesteps_all[[i]],
                  ".png"
                ),
                scale = 1.3
              )

              final_distribution_expanded_map_continuous_no_points_by_timestep
            }

  # rename time-specific list elements for clarity
  names(final_distribution_expanded_map_continuous_no_points_list) <-
    paste0("final_distribution_expanded_map_continuous_no_points_",
           timesteps_all)
}

{
  final_distribution_expanded_map_continuous_list <-
    foreach(i = seq_along(timesteps),
            .options.future = list(seed = TRUE)) %dofuture% {
              final_distribution_expanded_map_continuous_by_timestep <-
                {
                  boundaries_plot_layer + final_distribution_expanded_map_continuous_layer_list[[i]]
                }

              # save the image output
              tmap_save(
                final_distribution_expanded_map_continuous_by_timestep,
                filename = str_c(
                  "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/final_distribution_expanded_map_continuous_",
                  timesteps[[i]],
                  ".png"
                ),
                scale = 1.3
              )

              final_distribution_expanded_map_continuous_by_timestep
            }

  # rename time-specific list elements for clarity
  names(final_distribution_expanded_map_continuous_list) <-
    paste0("final_distribution_expanded_map_continuous_",
           timesteps)
}

{
  plan(sequential)
}  # end parallel processing


# Save all objects in R environment ----

# save all objects in R environment
save.image(
  "outputs/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil/C_variipennis_dynamic_nestedcv_reduced_vars_NAfilledSoil.RData"
)
