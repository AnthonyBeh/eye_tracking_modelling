# label and position information from text images
#
# ds 2020-08-18
library(dplyr)
library(tibble)
library(imager)
# BiocManager::install("EBImage")
library(EBImage)
library(tidyverse)
library(ggplot2)
library(imager)

# load image --------------------------------------------------------------
fname <- "E:\\thumbdrive\\gazeFilter\\textImg\\4.png"
textImage <- EBImage::readImage(fname)

# and make sure it's grayscale
textImage <- channel(textImage, "gray")

# display
EBImage::display(textImage)

# threshold
textImage = thresh(textImage, 10, 10, 0.05) 
EBImage::display(textImage)

# fill / erode / turn to blobs --------------------------------------------
textImage_blobbed <- closing(textImage, makeBrush(7, shape='disc'))
textImage_lablelled <-  bwlabel(textImage_blobbed)

EBImage::display(textImage_blobbed)
EBImage::display(textImage_lablelled)

# compute + convert to dataframe
fts <-  computeFeatures.shape(textImage_lablelled ) %>% 
  as_tibble()

# compute + convert to dataframe
fts <-  computeFeatures.shape(textImage_blobbed) %>% 
  as_tibble()

fts
# some useful metrics
glimpse(fts)

# the moments are central tendency of x, y of each labelled blob
# major axis, angle, etc. (assuming bivariate gaussian, I assume - check doc'n)
ftm <-  computeFeatures.moment(textImage_lablelled ) %>% 
  as_tibble()

glimpse(ftm)

# now plot moment data on image -------------------------------------------

d <- dim(textImage)

ftm %>% 
  ggplot(aes(x = m.cx, 
             y = m.cy, 
             angle = m.theta, 
             radius=m.majoraxis)) +
  geom_spoke() + # maybe not? not quite correct /// but to see what you can do
  geom_point(size=2, color="red", alpha=0.5) +
  scale_y_reverse(limits=c(d[2],0)) +
  scale_x_continuous(limits=c(0, d[1])) +
  labs(x="pixels in x", y = "pixels in y", title = "centroids and major axis / theta") +
  theme_minimal()

