##
##
## SCRIPT NAME: vsv_all_occurrence_mean_yearly_interval_by_state
##
##

### Set working environment ----

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

p_load(tidyverse,
       raster,
       terra,
       sf,
       tmap,
       ragg,
       viridis,
       pryr)

# get citation and version information about loaded packages
sink(
  "outputs/vsv_all_occurrence_mean_yearly_interval_by_state/loaded_packages_info.txt"
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
if (dir.exists("outputs/vsv_all_occurrence_mean_yearly_interval_by_state/raster_tempdir") == FALSE) {
  dir.create("outputs/vsv_all_occurrence_mean_yearly_interval_by_state/raster_tempdir")
}
raster_tempdir <-
  "outputs/vsv_all_occurrence_mean_yearly_interval_by_state/raster_tempdir/"


### Prepare presence locations ----

# import and check case occurrence data
occurrence <-
  read.csv(
    "data/occurrence_data/VSV_Mexico_merged_012621_fixed.csv",
    fileEncoding = readr::guess_encoding("data/occurrence_data/VSV_Mexico_merged_012621_fixed.csv")[[1, 1]]
  ) %>%
  filter(!is.na(Year)) %>%  # exclude records with no year
  filter(!is.na(Occurrence)) %>%  # exclude records with no confirmed occurrence of VS
  rename(
    year = Year,
    infection_type = Inf_.type,
    latitude = Lat,
    longitude = Long
  ) %>%
  mutate(State = str_to_title(State))
head(occurrence)
str(occurrence)
dim(occurrence)

# save data frame as csv file
write.csv(
  occurrence,
  file = "outputs/vsv_all_occurrence_mean_yearly_interval_by_state/occurrence.csv",
  row.names = FALSE,
  fileEncoding = readr::guess_encoding("data/occurrence_data/VSV_Mexico_merged_012621_fixed.csv")[[1, 1]]
)

### Calculate state mean yearly case intervals ----

# for each state, calculate the yearly interval between records

occurrence_year_diff <-
  occurrence %>% dplyr::select(State, year) %>% dplyr::arrange(State, year)
occurrence_year_diff

occurrence_year_diff$year_diff <-
  with(occurrence_year_diff,
       ave(
         as.numeric(year),
         State,
         FUN = function(x)
           c(0, diff(x))
       ))
occurrence_year_diff
summary(occurrence_year_diff)

# exclude first-year records for each state (since there are no earlier records for comparison to calculate their corresponding yearly interval)
occurrence_year_diff_first_year_removed <-
  occurrence_year_diff %>%
  group_by(State) %>%
  dplyr::filter(year != min(year))
occurrence_year_diff_first_year_removed

# for each state, change all values of zero to match the non-zero value per year (to correct the stepwise interval calculations)
occurrence_year_diff_first_year_removed_with_fixed_values <-
  occurrence_year_diff_first_year_removed %>%
  group_by(State, year) %>%
  dplyr::mutate(year_diff = max(year_diff))
occurrence_year_diff_first_year_removed_with_fixed_values

occurrence_year_diff_old <- occurrence_year_diff

occurrence_year_diff <-
  occurrence_year_diff_first_year_removed_with_fixed_values

# save data frame as csv file
write.csv(
  occurrence_year_diff,
  file = "outputs/vsv_all_occurrence_mean_yearly_interval_by_state/occurrence_year_diff.csv",
  row.names = FALSE,
  fileEncoding = readr::guess_encoding("data/occurrence_data/VSV_Mexico_merged_012621_fixed.csv")[[1, 1]]
)

# calculate mean yearly intervals for all multi-year states
occurrence_state_mean_yearly_intervals <- occurrence_year_diff %>%
  group_by(State) %>%
  summarize(mean_year_diff = mean(year_diff))
occurrence_state_mean_yearly_intervals
summary(occurrence_state_mean_yearly_intervals)

# isolate and summarize records from single-year states
occurrence_single_year_states <-
  occurrence_year_diff_old %>%
  group_by(State) %>%
  distinct(year, .keep_all = TRUE) %>%
  add_tally(name = "num_years") %>%
  dplyr::filter(num_years == 1) %>%
  summarize(mean_year_diff = mean(year_diff)) %>%
  dplyr::mutate(mean_year_diff = NA)

# insert single-year state rows into mean yearly intervals data frame
occurrence_state_mean_yearly_intervals_with_single_year_states <-
  occurrence_state_mean_yearly_intervals %>%
  dplyr::union(occurrence_single_year_states) %>%
  dplyr::arrange(State) #%>%
# dplyr::mutate(mean_year_diff = ifelse(mean_year_diff == 1000, "NA (only one year available)", mean_year_diff))
occurrence_state_mean_yearly_intervals_with_single_year_states

occurrence_state_mean_yearly_intervals_old <-
  occurrence_state_mean_yearly_intervals

occurrence_state_mean_yearly_intervals <-
  occurrence_state_mean_yearly_intervals_with_single_year_states
occurrence_state_mean_yearly_intervals

# save data frame as csv file
write.csv(
  occurrence_state_mean_yearly_intervals,
  file = "outputs/vsv_all_occurrence_mean_yearly_interval_by_state/occurrence_state_mean_yearly_intervals.csv",
  row.names = FALSE,
  fileEncoding = readr::guess_encoding("data/occurrence_data/VSV_Mexico_merged_012621_fixed.csv")[[1, 1]]
)


### Define study area ----

# load elevation raster
elevation <-
  raster::raster(
    str_c(
      main_path,
      "environmental_data/NEW_environmental_data/EarthEnv_Topography/final/elevation_4KMmd_GMTEDmd.tif"
    )
  )
elevation

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
    "data/study_areas/Ecoregions_Study_Area_VSV_Mexico_Endemic_All_Extended_Convex_Hull_Mask.gpkg"
  ) %>%
  st_transform(st_crs(elevation)) %>%
  st_make_valid()


### Mapping endemicity ----

# import SWD object (which includes filtered presence data) from main script as data frame
data <-
  read.csv("outputs/vsv_all_dynamic_nestedcv_final_reduced_vars/data.csv")

# define presence localities that were included in the SWD object (i.e., excluding points dropped due to missing environmental data)
presence_final <-
  data %>% dplyr::filter(pa == 1) %>%
  dplyr::select(X, Y) %>%
  rename(longitude = X,
         latitude = Y)

# make sf object versions for final presence coords
presence_final_sf <- st_as_sf(
  x = presence_final,
  coords = c("longitude", "latitude"),
  crs = st_crs(elevation)
)

boundaries_mexico_states <- read_sf(
  str_c(
    main_path,
    "environmental_data/study_areas/reference_boundaries/Global_States_Provinces.gpkg"
  )
) %>%
  st_transform(st_crs(elevation)) %>%
  st_make_valid() %>%
  filter(geonunit == "Mexico") %>%
  dplyr::select(name)

# add mean yearly interval data to mexico state boundaries object
boundaries_mexico_states_occ_mean_yearly_intervals <-
  boundaries_mexico_states %>%
  full_join(occurrence_state_mean_yearly_intervals, by = c("name" = "State"))

# create plot layer from country and state/province boundaries
boundaries_alt_plot_layer <- {
  tm_shape(
    boundaries_countries,
    projection = st_crs(boundaries_mexico_states_occ_mean_yearly_intervals),
    bbox = st_bbox(boundaries_mexico_states_occ_mean_yearly_intervals)
  ) +
    tm_polygons(alpha = 0, border.alpha = 0) +  # transparent copy of polygon (baseline)
    tm_graticules(alpha = 0.2,
                  ticks = FALSE) +
    tm_layout(inner.margins = 0,
              bg.color = "white") +  # add graticules (long-lat lines) first
    tm_shape(
      boundaries_countries,
      projection = st_crs(boundaries_mexico_states_occ_mean_yearly_intervals),
      bbox = st_bbox(boundaries_mexico_states_occ_mean_yearly_intervals)
    ) +
    tm_polygons(col = "white", border.alpha = 0) +  # opaque grey copy of polygon (to hide graticules under polygon)
    tm_shape(
      boundaries_countries,
      projection = st_crs(boundaries_mexico_states_occ_mean_yearly_intervals),
      bbox = st_bbox(boundaries_mexico_states_occ_mean_yearly_intervals)
    ) +
    tm_polygons(
      alpha = 0.15,
      border.alpha = 1,
      border.col = "black",
      lwd = 1
    ) +  # visible copy of country boundaries
    tm_shape(
      boundaries_states_provinces,
      projection = st_crs(boundaries_mexico_states_occ_mean_yearly_intervals),
      bbox = st_bbox(boundaries_mexico_states_occ_mean_yearly_intervals)
    ) +
    tm_polygons(
      alpha = 0,
      border.alpha = 1,
      border.col = "black",
      lwd = 0.8
    )  # visible copy of state/province boundaries
}

# create mean yearly interval maps

mexico_states_occ_mean_yearly_intervals_map_no_points_layer <-
  {
    tm_shape(
      boundaries_mexico_states_occ_mean_yearly_intervals,
      projection = st_crs(boundaries_mexico_states_occ_mean_yearly_intervals),
      bbox = st_bbox(boundaries_mexico_states_occ_mean_yearly_intervals)
    ) + tm_polygons(
      col = "mean_year_diff",
      palette = viridis(
        256,
        option = "magma",
        begin = 0.05,
        end = 0.95,
        direction = -1
      ),
      alpha = 1,
      title = "Mean yearly interval\nbetween VS cases\n(1981–2020)",
      legend.reverse = TRUE,
      interval.closure = "left",
      legend.format = list(scientific = TRUE, format = "f"),

      # # VARIANT 1
      # style = "fixed",
      # breaks = seq(1, max(
      #   summary(
      #     boundaries_mexico_states_occ_mean_yearly_intervals$mean_year_diff
      #   )
      # ), 1)

      # # VARIANT 2
      # style = "fixed",
      # breaks = seq(1, max(
      #   summary(
      #     boundaries_mexico_states_occ_mean_yearly_intervals$mean_year_diff
      #   )
      # ), 2)

      # VARIANT 3
      style = "fixed",
      breaks = c(seq(1, 2, 1),
                 seq(3, max(
                   summary(
                     boundaries_mexico_states_occ_mean_yearly_intervals$mean_year_diff
                   )
                 ), 2))

      # # VARIANT 4
      # style = "cont"

      # # VARIANT 5
      # style = "cont",
      # breaks = seq(1, max(
      #   summary(
      #     boundaries_mexico_states_occ_mean_yearly_intervals$mean_year_diff
      #   )
      # ), 2)

    ) + tm_legend(legend.outside = TRUE)
  }
mexico_states_occ_mean_yearly_intervals_map_no_points_layer

mexico_states_occ_mean_yearly_intervals_map_layer <- {
  mexico_states_occ_mean_yearly_intervals_map_no_points_layer + tm_shape(presence_final_sf) + tm_dots(
    size = 0.08,
    shape = 21,
    col = "blue",
    alpha = 0.7,
    border.lwd = 0.55,
    border.alpha = 0.9
  )
}

(
  mexico_states_occ_mean_yearly_intervals_map <-
    boundaries_alt_plot_layer + mexico_states_occ_mean_yearly_intervals_map_layer
)

(
  mexico_states_occ_mean_yearly_intervals_map_no_points <-
    boundaries_alt_plot_layer + mexico_states_occ_mean_yearly_intervals_map_no_points_layer
)

# save the image output
tmap_save(
  mexico_states_occ_mean_yearly_intervals_map,
  filename = "outputs/vsv_all_occurrence_mean_yearly_interval_by_state/mexico_states_occ_mean_yearly_intervals_map_v3.png",
  scale = 1.7
)

# save the image output
tmap_save(
  mexico_states_occ_mean_yearly_intervals_map_no_points,
  filename = "outputs/vsv_all_occurrence_mean_yearly_interval_by_state/mexico_states_occ_mean_yearly_intervals_map_no_points_v3.png",
  scale = 1.7
)


#### With study area overlay ----

# build an larger extent for next set of maps
bbox_extended <-
  st_bbox(boundaries_mexico_states_occ_mean_yearly_intervals)

bbox_extended[4] <- st_bbox(study_area_mask)[4]

# create plot layer from country and state/province boundaries
boundaries_alt_with_study_area_plot_layer <- {
  tm_shape(
    boundaries_countries,
    projection = st_crs(bbox_extended),
    bbox = st_bbox(bbox_extended)
  ) +
    tm_polygons(alpha = 0, border.alpha = 0) +  # transparent copy of polygon (baseline)
    tm_graticules(alpha = 0.2,
                  ticks = FALSE) +
    tm_layout(inner.margins = 0,
              bg.color = "white") +  # add graticules (long-lat lines) first
    tm_shape(
      boundaries_countries,
      projection = st_crs(bbox_extended),
      bbox = st_bbox(bbox_extended)
    ) +
    tm_polygons(col = "white", border.alpha = 0) +  # opaque grey copy of polygon (to hide graticules under polygon)
    tm_shape(
      boundaries_countries,
      projection = st_crs(bbox_extended),
      bbox = st_bbox(bbox_extended)
    ) +
    tm_polygons(
      alpha = 0.15,
      border.alpha = 1,
      border.col = "black",
      lwd = 1
    ) +  # visible copy of country boundaries
    tm_shape(
      boundaries_states_provinces,
      projection = st_crs(bbox_extended),
      bbox = st_bbox(bbox_extended)
    ) +
    tm_polygons(
      alpha = 0,
      border.alpha = 1,
      border.col = "black",
      lwd = 0.8
    )  # visible copy of state/province boundaries
}

# create plot layer from study area mask
study_area_mask_alt_plot_layer <- {
  tm_shape(
    study_area_mask,
    projection = st_crs(bbox_extended),
    bbox = st_bbox(bbox_extended),
    is.master = TRUE
  ) +
    tm_polygons(
      col = "cornsilk",
      alpha = 0.2,
      border.alpha = 0.5,
      lwd = 1.8
    )
}

# create mean yearly interval maps

(
  mexico_states_occ_mean_yearly_intervals_map_with_study_area <-
    boundaries_alt_with_study_area_plot_layer + mexico_states_occ_mean_yearly_intervals_map_layer + study_area_mask_alt_plot_layer
)

(
  mexico_states_occ_mean_yearly_intervals_map_with_study_area_no_points <-
    boundaries_alt_with_study_area_plot_layer + mexico_states_occ_mean_yearly_intervals_map_no_points_layer +
    study_area_mask_alt_plot_layer
)

# save the image output
tmap_save(
  mexico_states_occ_mean_yearly_intervals_map_with_study_area,
  filename = "outputs/vsv_all_occurrence_mean_yearly_interval_by_state/mexico_states_occ_mean_yearly_intervals_map_with_study_area_v3.png",
  scale = 1.7
)

# save the image output
tmap_save(
  mexico_states_occ_mean_yearly_intervals_map_with_study_area_no_points,
  filename = "outputs/vsv_all_occurrence_mean_yearly_interval_by_state/mexico_states_occ_mean_yearly_intervals_map_with_study_area_no_points_v3.png",
  scale = 1.7
)


## Save all objects in R environment ----

# save all objects in R environment
save.image(
  "outputs/vsv_all_occurrence_mean_yearly_interval_by_state/vsv_all_occurrence_mean_yearly_interval_by_state.RData"
)

