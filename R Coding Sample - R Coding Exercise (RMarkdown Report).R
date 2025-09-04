#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------
#' ***Melanie Veron***
#' 
#' 
#' ***October 19, 2017***
#' 
#' 
#' ***Introduction to R 3: Working With Data***
#' 
#' 
#' ***In-Class Exercise***
#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------

library(ggplot2)
library(dplyr)

# clear history
rm(list = ls())

# import birdsticks2007 dataset
# assign to tick.dat
tick.dat <- read.csv("C:\\Users\\Melanie\\Documents\\School\\Postgrad\\Glasgow\\Academics\\Courses\\Key Research Skills\\Lecture Notes - Intro to R\\Intro to R 3 - Working With Data\\birdsticks2007.csv")

# check dataset
glimpse(tick.dat)
head(tick.dat)

# see number of individuals per species trapped each month
table(tick.dat$species, tick.dat$month)

#-----------------------------------------------------------------------------------
#' **1. Adding and Redefining Variables**
#-----------------------------------------------------------------------------------

#' a) Adding a new variable

# wt.bag: weight of bird included in handling bag
# bagwt: weight of bag (if used)
# bird.wt: actual weight of bird
#
# make new bird.wt variable
# by subtracting bagwt from wt.bag
tick.dat <- mutate(tick.dat, bird.wt = wt.bag - bagwt)

# check that dataset has bird.wt variable
glimpse(tick.dat)

#' b) Recording an existing variable

# retrap: shows whether individual has been caught before ('1') 
# or caught for first time ('0')
#
# use recode() to redefine retrap as a factor variable 
# and show two factor levels as 'yes' (replacing '1')
# and 'no' (replacing '0')
# i.e., replace numerical retrap variable with categorical retrap variable
#
# must do this using mutate()
# first define dataset (tick.dat)
# then transforming retrap via recode()
# before assigning result to tick.dat
tick.dat <- mutate(tick.dat, retrap = recode(retrap, '1' = "yes", '0' = "no"))

# check that dataset has retrap as a factor variable
glimpse(tick.dat)

#' c) Creating a new data variable

# make new date variable with day.month.year format
# to show full date for each observation
# rather than having just three separate columns 
# for day, month, and year per observation
#
# must do this using mutate()
# first define dataset (tick.dat)
# then with paste(), copy values from day, month, and year variables
# and paste them into a single new date variable
# with a '.' separating each value per slot
# then assign result to tick.dat 
tick.dat <- mutate(tick.dat, date = paste(day, month, year, sep = "."))

# check dataset for date variable
glimpse(tick.dat)

# convert date variable into dates recognized by R
# by using dmy() from lubridate package
library(lubridate)

# dmy(): converts values to dates that R recognizes
# (in year-month-day format, e.g., 2007-01-27)
tick.dat <- mutate(tick.dat, date = dmy(date)) 

# check that date variable is now R-friendly
glimpse(tick.dat)

#-----------------------------------------------------------------------------------
#' **2. Sorting**
#-----------------------------------------------------------------------------------

# using arrange(), sort data by species then by sex in ascending order 
tick.dat <- arrange(tick.dat, species, sex)

# check dataset for new order
select(tick.dat, species, sex)

# same, but now in descending order
tick.dat <- arrange(tick.dat, desc(species, sex))

# check dataset for new order
select(tick.dat, species, sex)

# sort data by chronological order (i.e., by date)
tick.dat <- arrange(tick.dat, date)

# check dataset for new order
select(tick.dat, date)

#-----------------------------------------------------------------------------------
#' **3. Subsetting Your Data**
#-----------------------------------------------------------------------------------

# make subset of tick.dat
# including observations only from April though October
# by using filter()
# and setting filter to be 4 <= month <= 10 from tick.dat
# then assign result to new data frame tick.season
tick.season <- filter(tick.dat, month >= 4 & month <= 10)

# check tick.season dataset
glimpse(tick.season)

#-----------------------------------------------------------------------------------
#' **4. Dealing With NAs**
#-----------------------------------------------------------------------------------

# calculate median number of ticks (ticktotal) for each species
# 
# via piping method, use group_by() to group by species within tick.season
# then use summarise() to determine median of ticktotal
# and assign those medians per species to median.tick variable
#
# notice in the results how some of the species have missing medians
# due to missing ticktotal values (NAs) from those species groups

tick.season %>%
  group_by(species) %>%
      summarise(median.tick = median(ticktotal)) 

# exclude rows where there are missing values (NAs) 
# in ticktotal variable from tick.season data frame
#
# use is.na() to select all NAs in ticktotal from tick.season
# and use '!' action before that to exclude all results 
# (that are to right-hand side of '!')
# from tick.season 
# then assign result to tick.season
tick.season <- tick.season[!is.na(tick.season$ticktotal),] 

# overwrite default where printing data frame results
# shows only first ten rows
options(tibble.print_max = Inf)

# check that there are no NAs in ticktotal
filter(tick.season, is.na(tick.season$ticktotal))

# again, calculate median ticktotal for each species
tick.season %>%
  group_by(species) %>%
    summarise(median.tick = median(ticktotal))

# assign median ticktotal result to new object tick.summary
tick.summary <- tick.season %>%
                  group_by(species) %>%
                   summarise(median.tick = median(ticktotal))

# sort tick.summary data frame by median.tick in descending order
# (so that species with the most ticks come first)
tick.summary <- arrange(tick.summary, desc(median.tick))

# check tick.summary dataset
glimpse(tick.summary)

#-----------------------------------------------------------------------------------
#' **5. Plotting the Data**
#-----------------------------------------------------------------------------------

# make point plot with species (x-axis) vs. median.tick (y-axis)
# and change x-axis labels (species names) so that they're readable
# 
# theme() line changes the x-axis text from being parallel to the x-axis
# (which makes the labels unreadable 
# since the species names side-by-side run into each other)
# to being perpendicular (more legible)
ggplot(tick.summary, aes(x = species, y = median.tick)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90))

# the plot ends up not being sorted by median.tick
# since ggplot by default wants the factor levels
# in a categorical variable (i.e., species) in alphabetical order
#
# to fix this, change species variable so that it is ordered by median.tick
# rather than by the default alphabetical order
tick.summary$species <- factor(tick.summary$species, levels = 
                          tick.summary$species[rev(order(tick.summary$median.tick))])

# plot the same graph and check that species are sorted by median.tick
ggplot(tick.summary, aes(x = species, y = median.tick)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90))

# plot species (x-axis) vs. ticktotal (y-axis) from tick.season
# use boxplot to see not just the median tick numbers
# but also how the data are distributed
ggplot(tick.season, aes(x = species, y = ticktotal)) +
  geom_boxplot() +
  geom_point() +
  xlab("Bird species") +
  ylab("Number of ticks") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))


