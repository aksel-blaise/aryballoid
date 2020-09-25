# load packages
library(here)
library(Momocs)

# read images
jpg.list <- list.files(here("/img.perdiz"), full.names = TRUE)

# read attribute data
att.data <- read.csv("perdiz.csv", header = TRUE, as.is = TRUE)

# attribute to factor
att.data$type <- as.factor(att.data$type)
att.data$trinomial <- as.factor(att.data$trinomial)
att.data$raw.mat <- as.factor(att.data$raw.mat)

# generate outlines
outlines <- jpg.list %>%
  import_jpg()

# add attributes
data.out <- Out(outlines, 
                fac = att.data)

# scale, align, and center specimens
norm.outlines <- data.out %>% 
  coo_scale() %>%
  coo_align() %>% 
  coo_center()

# outline pile
pile(norm.outlines)

panel(norm.outlines, names = TRUE)
# mosaic of individual specimens from the different sites
mosaic(norm.outlines, ~trinomial)
# mosaic of individual specimens rendered from different materials
mosaic(norm.outlines, ~raw.mat)

# calibrate how many harmonics needed
calibrate_harmonicpower_efourier(norm.outlines, 
                                 nb.h = 30)

# 11 harmonics needed to capture 99 percent of variation
calibrate_reconstructions_efourier(norm.outlines, 
                                   range = 1:11)

# generate efa outlines with 11 harmonics
efa.outlines <- efourier(norm.outlines, 
                         nb.h = 11, 
                         norm = TRUE)

# use efa.outlines for pca
pca.outlines <- PCA(efa.outlines)

# pca 
scree_plot(pca.outlines)

# plot pca by site
plot_PCA(pca.outlines, 
         morphospace_position = "range",
         ~trinomial, zoom = 1.3)

# plot pca by raw material
plot_PCA(pca.outlines, 
         morphospace_position = "range",
         ~raw.mat, zoom = 1.3)

# contribution of each pc
# by site
boxplot(pca.outlines, ~trinomial, nax = 1:5)
# by raw material
boxplot(pca.outlines, ~raw.mat, nax = 1:5)

# mean shape + 2sd for the first 5 pcs
PCcontrib(pca.outlines, nax = 1:5)

# manova

# shape difference between sites?
MANOVA(pca.outlines, 'trinomial')
# which differ?
MANOVA_PW(pca.outlines, 'trinomial')

# shape difference between raw.mat?
MANOVA(pca.outlines, 'raw.mat')
# which differ?
MANOVA_PW(pca.outlines, 'raw.mat')

# mean shapes

# site
ms.2 <- MSHAPES(efa.outlines, ~trinomial)
plot_MSHAPES(ms.2, size = 0.8)

# raw material
ms.3 <- MSHAPES(efa.outlines, ~raw.mat)
plot_MSHAPES(ms.3, size = 0.75)

#end of code
