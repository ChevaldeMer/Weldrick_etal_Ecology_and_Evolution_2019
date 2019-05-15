# SCRIPT 1

## Test SIA values for normality

## Read in data
SIAmanova <- read.csv("MANOVA.csv", header = TRUE)

### Shapiro-Wilks test for normality
shapiro.test(SIAmanova$d13C)
shapiro.test(SIAmanova$d15N)

### Levene's test for homogeneity of variance
library(car)
leveneTest(d13C ~ species, data = SIAmanova)
leveneTest(d15N ~ species, data = SIAmanova)

## MANOVA
d13C <- SIAmanova$d13C
d15N <- SIAmanova$d15N

res.man <- manova(cbind(d13C, d15N) ~ lat + lon + depth + species, 
                  data = SIAmanova)
summary(res.man, test = "Wilks")

Y <- cbind(d13C, d15N)
lat <- SIAmanova$lat
lon <- SIAmanova$lon
depth <- SIAmanova$depth
species <- SIAmanova$species
fit <- manova(Y ~ lat*lon*depth*species)
summary(fit, test = "Wilks")
summary(fit, test = "Pillai")
summary(fit, test = "Hotelling-Lawley")
summary(fit, test = "Roy")
summary.aov(fit, test = "Wilks")
summary.aov(fit, test = "Pillai")
summary.aov(fit, test = "Hotelling-Lawley")
summary.aov(fit, test = "Roy")

## Figure 2
## Abundance map (Fig 2) for relative pteropod species and POM

library(marmap)
library(ggplot2)
library(rgdal)
library(ggmap)
library(GEOmap)
library(gridExtra)
library(gtable)
library(grid)
library(rworldmap)
library(dplyr)
library(geosphere)
library(gpclib)
library(ggthemes)


## Obtain bathymetric data for map
kaxis <- getNOAA.bathy(lon1 = 60, lon2 = 100,
                       lat1 = -70, lat2 = -50, resolution = 1, keep = TRUE)

## Read in abundance data
setwd("~/Desktop/OneDrive - University of Tasmania/PhD/Chapters/Project 2 - Pteropod SIA")
### pteropods
abundkaxis <- read.csv("abundance bubble map.csv", header = TRUE)
### POM
POM <- read.csv("POM SIA results.csv", header = TRUE)

## Read in Bestley et al. 2018 frontal zone data
eddy <- read.csv("KAXIS_FEATURES_EDDY.csv", header = T)
saccf <- read.csv("KAXIS_FEATURES_SACCF.csv", header = T)
saf <- read.csv("KAXIS_FEATURES_ASF.csv", header = T)
sbdy <- read.csv("KAXIS_FEATURES_SBDY.csv", header = T)
ftj <- read.csv("KAXIS_FEATURES_FTJ.csv", header = T)
winter <- read.csv("KAXIS_FEATURES_WINTER.csv", header = T)
summer <- read.csv("KAXIS_FEATURES_SUMMER.csv", header = T)

## Read in coastal features
coast <- readOGR(getwd(),"all_coast_poly_2003")
coast_df <- fortify(coast)
icebergs <- readOGR(getwd(),"icebergs_dataset_309")
icebergs_df <- fortify(icebergs)
tongue <- readOGR(getwd(), "ice_tongues_boundaries_dataset_309")
tongue_df <- fortify(tongue)
ice <- readOGR(getwd(), "ice_poly_2003")
ice_df <- fortify(ice)


## Customise colours for each pteropod species
cols <- c("Clio pyramidata sulcata" = "#CCBA71",
          "Clione limacina antarctica" = "#0E0D0D",
          "Gymnosome spp." = "#00A08A",
          "Spongiobranchaea australis" = "#9985A5")

## Customise labels for taxa
mylabels <- list(expression(italic("Clio pyramidata sulcata")),
                 expression(italic("Clione limacina antarctica")),
                 "Gymnosome spp.",
                 expression(italic("Spongiobranchaea australis")))

## Functions for labelling latitude and longitude
scale_x_longitude <- function(xmin = -180, xmax = 180, step = 1, ...) {
  xbreaks <- seq(xmin,xmax,step)
  xlabels <- unlist(lapply(xbreaks, function(x) ifelse(x < 0, parse(text = paste0(abs(dms(x)$d),"^o", "*W")), ifelse(x > 0, parse(text = paste0(abs(dms(x)$d),"^o", "*E")),abs(dms(x))))))
  return(scale_x_continuous("Longitude", breaks = xbreaks, labels = xlabels, expand = c(0, 0), ...))
}
scale_y_latitude <- function(ymin = -90, ymax = 90, step = 0.5, ...) {
  ybreaks <- seq(ymin,ymax,step)
  ylabels <- unlist(lapply(ybreaks, function(x) ifelse(x < 0, parse(text=paste0(abs(dms(x)$d),"^o", "*S")), ifelse(x > 0, parse(text = paste0(abs(dms(x)$d),"^o", "*N")),abs(dms(x))))))
  return(scale_y_continuous("Latitude", breaks = ybreaks, labels = ylabels, expand = c(0, 0), ...))
}

## Main map
main <- autoplot(kaxis, geom = c("raster", "contour"), 
               colour = "white", size=0.1) +
  scale_fill_gradient2(low = "lightsteelblue3", 
                       mid = "mintcream", high = "whitesmoke", 
                       guide = FALSE) +
  geom_path(aes(x = lon, y = lat), data = eddy, linetype = 2, alpha = 0.4) +
  geom_path(aes(x = lon, y = lat), data = saccf, linetype = 2, alpha = 0.4) +
  geom_path(aes(x = lon, y = lat), data = saf, linetype = 2, alpha = 0.4) +
  geom_path(aes(x = lon, y = lat), data = sbdy, linetype = 2, alpha = 0.4) +
  geom_path(aes(x = lon, y = lat), data = ftj, linetype = 2, alpha = 0.4) +
  geom_path(aes(x = lon, y = lat), data = winter, linetype = 3, alpha = 0.4) +
  geom_path(aes(x = lon, y = lat), data = summer, linetype = 3, alpha = 0.4) +
  geom_path(aes(x = lon, y = lat), data = wfile, colour = "gray30") +
  geom_point(aes(x = longitude, y = latitude, 
                 size = abundance, colour = species), 
             data = abundkaxis, alpha = 0.8) +
  geom_point(aes(x = longitude, y = latitude, shape = fraction),
             data = POM) +
  scale_color_manual("species",
                     values = cols,
                     breaks = c("Clio pyramidata sulcata", 
                                "Clione limacina antarctica",
                                "Gymnosome spp.",
                                "Spongiobranchaea australis"),
                     labels = mylabels) +
  scale_size_continuous(breaks = c(10,50,100,150,250,270),
                        limits = c(1,300),
                        labels = c("<10","11-50","51-100","101-150","151-250",">250")) +
  scale_shape_manual(values = c(17,18)) +
  geom_polygon(data = coast_df,
               aes(x = long, y = lat, group = group),
               color = "snow4", fill = "whitesmoke", size = 0.25) +
  geom_polygon(data = icebergs_df,
               aes(x = long, y = lat, group = group),
               color = "snow4", fill = "azure", size = 0.25) +
  geom_polygon(data = tongue_df,
               aes(x = long, y = lat, group = group),
               color = "snow4", fill = "azure", size = 0.25) +
  geom_polygon(data = ice_df,
               aes(x = long, y = lat, group = group),
               color = "snow4", fill = "snow", size = 0.25) +
  annotate("text", x = 90, y = -68, label = c("Antarctica"),
           size = 3, colour = "gray18", fontface = "bold") +
  annotate("text", x = 64, y = -62.2, label = c("sACCf"), angle = -15, size = 3, colour = "gray44") +
  annotate("text", x = 64, y = -63.7, label = c("SB"), size = 3, colour = "gray44") +
  annotate("text", x = 82, y = -55.9, label = c("FTC"), size = 3, colour = "gray44") +
  annotate("text", x = 81.5, y = -65.14, label = c("ASF"), size = 3, colour = "gray44") +
  annotate("text", x = 76.5, y = -57.5, label = c("Kerguelen \n Plateau"), colour = "gray42") +
  labs(aes(x = "Longitude", y = "Latitude"), data = abundkaxis) +
  guides(size = guide_legend(nrow = 1, byrow=T), 
         color = guide_legend(nrow = 2, byrow = T),
         shape = guide_legend(title = "POM fraction", nrow = 1, byrow = T)) +
  coord_cartesian(xlim = c(62, 98), ylim = c(-69,-55)) +
  scale_x_longitude(60, 100, 5) +
  scale_y_latitude(-70, -55, 5) +
  theme_bw() +
  theme(panel.border = element_rect(fill = NA, colour = "black", size = 1),
        axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        legend.position = "bottom", 
        legend.title = element_blank(), 
        legend.box = "vertical",
        legend.text = element_text(size = 7), 
        legend.spacing.y = unit(-0.2, "cm")) 
main

## Set up inset map (adjusted from http://egallic.fr/maps-with-r/)
xlim <- range(abundkaxis$longitude) + c(-1, 0.5)
ylim <- range(abundkaxis$latitude) + c(-0.5, 1)

## World map
worldMap <- getMap()
world.points <- fortify(worldMap)
world.points$region <- world.points$id

world.df <- world.points[,c("long","lat","group", "region")]

worldmap <- ggplot() + 
  geom_polygon(data = world.df, aes(x = long, y = lat, group = group)) +
  scale_y_continuous(breaks = (-2:2) * 30) +
  scale_x_continuous(breaks = (-4:4) * 45)

worldmap

## Change projection of inset map
worldmap <- ggplot() + 
  geom_polygon(data = world.df, aes(x = long, y = lat, group = group)) + theme_bw() +
  geom_rect(data = data.frame(),
            aes(xmin = xlim[1], xmax = xlim[2], ymin = ylim[1], ymax = ylim[2]),
            colour = "red", fill = NA) +
  scale_y_continuous(breaks = (-2:2) * 30) +
  scale_x_continuous(breaks = (-4:4) * 45) +
  coord_map("ortho", orientation=c(-90, 0, 0)) +
  labs(x = NULL, y = NULL) + theme_map() + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_line(colour = "gray65"))
worldmap

## Combine maps
grid.newpage()
v1 <- viewport(width = 1, height = 1, x = 0.5, y = 0.5)  # the larger map
v2 <- viewport(width = 0.25, height = 0.25, x = 0.87, y = 0.85)  # the inset map
print(main, vp = v1)
print(worldmap, vp = v2)

## Figure 3
## This code script for plotting the isotopic niche widths and densities of co-occurring pteropods and POM

## Read in isotopes data
### For pteropods
dragon <- read.csv("dragonplot.csv", header = TRUE)
### For POM
POM <- read.csv("POM SIA results.csv", header = TRUE)

library(ggplot2)
library(siar)
library(plyr)
library(gridExtra)
library(grid)
library(png)
library(RCurl)

## Subset data manually
sp1 <- subset(dragon, dragon$species=="Clio pyramidata f. sulcata")
sp2 <- subset(dragon, dragon$species=="Clione limacina antarctica")
sp3 <- subset(dragon, dragon$species=="Spongiobranchaea australis")
sp4 <- subset(dragon, dragon$species=="Large-fraction POM")
sp5 <- subset(dragon, dragon$species=="Small-fraction POM")

## Calculate standard.ellipse for each subset
SE1 <- standard.ellipse(sp1$d13C, sp1$d15N)
SE2 <- standard.ellipse(sp2$d13C, sp2$d15N)
SE3 <- standard.ellipse(sp3$d13C, sp3$d15N)
SE4 <- standard.ellipse(sp4$d13C, sp4$d15N)
SE5 <- standard.ellipse(sp5$d13C, sp5$d15N)

## Create name variable for each group (so ggplot2 knows how to colour it)
## The name needs to be the same as original dataframe containing isotopic data
sp1_ <- rep("Clio pyramidata f. sulcata", length(SE1$xSEAc))
sp2_ <- rep("Clione limacina antarctica", length(SE2$xSEAc))
sp3_ <- rep("Spongiobranchaea australis", length(SE3$xSEAc))
sp4_ <- rep("Large-fraction POM", length(SE4$xSEAc))
sp5_ <- rep("Small-fraction POM", length(SE5$xSEAc))

## Create new dataframe with names and ellipse outputs
species <- c(sp1_,sp2_,sp3_,sp4_,sp5_)
x <- c(SE1$xSEAc,SE2$xSEAc,SE3$xSEAc,SE4$xSEAc,SE5$xSEAc)
y <- c(SE1$ySEAc,SE2$ySEAc,SE3$ySEAc,SE4$ySEAc,SE5$ySEAc)
df_SE <- data.frame(x,y,species)
plot(df_SE$x, df_SE$y)

## Customise colours
cols <- c("Clio pyramidata f. sulcata" = "#CCBA71",
          "Clione limacina antarctica" = "#0E0D0D",
          "Spongiobranchaea australis" = "#9985A5",
          "Large-fraction POM" = "#F58231",
          "Small-fraction POM" = "royalblue4")

## Customise labels
mylabels <- list(expression(italic("C. pyramidata")),
                 expression(italic("C. antarctica")),
                 expression(italic("S. australis")),
                 "large-fraction POM",
                 "small-fraction POM")

## Calculate hulls
find_hull <- function(dragon) dragon[chull(dragon$d13C, dragon$d15N), ]
hulls <- ddply(dragon, "species", find_hull)

## Plot the scatter
scatterplot <- ggplot(dragon, aes(x = d13C, y = d15N, color = species)) +
  theme_bw() +
  geom_polygon(data = hulls, alpha = 0.05, linetype=3) +
  geom_point(size = 1) +
  scale_colour_manual("species",
                      values = cols,
                      breaks = c("Clio pyramidata f. sulcata",
                                 "Clione limacina antarctica",
                                 "Spongiobranchaea australis",
                                 "Large-fraction POM",
                                 "Small-fraction POM"),
                      labels = mylabels) +
  theme(legend.position = "none") +
  ylab(expression(paste(delta^{15}, "N (\u2030)"))) +
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) +
  geom_path(data = df_SE, aes(x = x, y = y, color = species), linetype = 1, size = 1)
scatterplot

## Read image icons
icon1 <- readPNG("Clio pyramidata_yellow.png")
icon2 <- readPNG("Spongiobranchaea australis_purple.png")
icon3 <- readPNG("Clione limacina antarctica_grey.png")

## Plot image icons
a1 <- annotation_raster(icon1, ymin = 2, ymax = 3.5, xmin = -25, xmax=-24) #clio
a2 <- annotation_raster(icon2, ymin = 5, ymax = 6.4, xmin = -23.6, xmax = -23) #spongio
a3 <- annotation_raster(icon3, ymin = 4.5, ymax = 6, xmin = -28.6, xmax = -27.7) #clione
scatterplot2 <- scatterplot + a1 + a2 + a3

## Plot marginal density of x (along top panel)
xdensity <- ggplot(dragon, aes(d13C, fill = species)) +
  geom_density(alpha = .95) +
  scale_fill_manual("species",
                    values = cols,
                    breaks = c("Clio pyramidata f. sulcata",
                               "Clione limacina antarctica",
                               "Spongiobranchaea australis",
                               "Large-fraction POM",
                               "Small-fraction POM"),
                    labels = mylabels) +
  theme(legend.position = "none",
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank())
xdensity

# Plot marginal density of y (along right panel)
ydensity <- ggplot(dragon, aes(d15N, fill = species)) +
  geom_density(alpha = .95) +
  scale_fill_manual("species",
                    values = cols,
                    breaks = c("Clio pyramidata f. sulcata",
                               "Clione limacina antarctica",
                               "Spongiobranchaea australis",
                               "Large-fraction POM",
                               "Small-fraction POM"),
                    labels = mylabels) +
  theme(legend.position = "none",
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank()) +
  coord_flip()
ydensity 

## Create Stand-alone legend
g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

coverage_plot <- ggplot(data=dragon, aes(x = d13C, y = d15N, group = species, color = species)) + 
  geom_line(size = 1) + 
  scale_color_manual("species",
                     values = cols,
                     breaks = c("Clio pyramidata f. sulcata",
                                "Clione limacina antarctica",
                                "Spongiobranchaea australis",
                                "Large-fraction POM",
                                "Small-fraction POM"),
                     labels = mylabels) +
  geom_point(aes(colour = species), show.legend = T, size=3) +
  scale_x_discrete(labels = seq(1, 30.0, by=1)) +
  theme_bw() +
  theme(legend.position = "right", axis.text.x = element_text(angle = 90),
        legend.title = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_blank()) +
  labs(x = "carbon") +
  scale_shape_discrete() +
  guides(shape=guide_legend(override.aes = list(size = 2, linetype = 0)))

mylegend <- g_legend(coverage_plot)
p3 <- grid.draw(mylegend)

## Arrange plots together
grid.arrange(xdensity, mylegend, scatterplot2, ydensity,
             ncol = 2, nrow = 2, widths = c(4, 1.4), heights = c(1.4,4))

## Figure 4
## This code script for plotting the body size vs trophic positions

## Read in body lengths and trophic position data
trophicbody <- read.csv("body length vs d15N.csv", header = TRUE)
mm <- trophicbody$Length..mm.

# Create biplot 
biplot <- ggplot(trophicbody, aes(x = mm, y = tp, shape = species, colour = species)) +
  geom_point(shape = 16, size = 4, alpha = .4) +
  scale_colour_manual(values = c("#0E0D0D","#CCBA71", "#9985A5")) +
  xlab("Length (mm)") + ylab("Trophic position") +
  geom_smooth(method = "lm", fill = NA) +
  theme_minimal() + theme(legend.position = "none") +
  theme(legend.justification = c(1,1), legend.position = c(1,1),
        legend.text = element_text(face = "italic"),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "transparent",
                                         colour = "transparent"),
        legend.box.background = element_rect(colour = "grey"))
biplot

## Read image icons
icon1 <- readPNG("Clio pyramidata_yellow.png")
icon2 <- readPNG("Clione limacina antarctica_grey.png")
icon3 <- readPNG("Spongiobranchaea australis_purple.png")

## Add icons to plot
a1 <- annotation_raster(icon1, ymin = 2.1, ymax = 2.5, xmin = 18.2, xmax=20) #clio
a2 <- annotation_raster(icon2, ymin = 2.75, ymax = 3.15, xmin = 18, xmax = 19.6) #clione
a3 <- annotation_raster(icon3, ymin = 3.3, ymax = 3.7, xmin = 9.8, xmax = 11) #spongio

biplot_icons <- biplot + a1 + a2 + a3
biplot_icons

## The following script feature results of linear regression analyses for each species
clio <- lm(d15N ~ mm, data = cliolengths)
summary(clio) #y = 4.74 - 0.08x + e where e ~ N(0, 0.60); R^2=0.17; p = 0.06

clionelengths <- read.csv("clione body lengths.csv", header = TRUE)
clione <- lm(d15N ~ mm, data = clionelengths)
summary(clione) #y = 6.06 - 0.11x + e where e ~ N(0, 0.70); R^2=0.18; p = 0.09

spongiolengths <- read.csv("spongio body lengths.csv", header = TRUE)
spongio <- lm(d15N ~ mm, data = spongiolengths)
summary(spongio) #y = 3.87 + 0.10x + e where e ~ N(0, 0.39); R^2= 0.39; p = 0.05

## Figure 5
## The following code creates the trophic position plots

library(tRophicPosition)
library(reshape2)
library(ggridges)
library(ggplot2)
library(grid)
library(png)
library(RCurl)

## Read in data
pterosTP <- read.csv("trophicposition pteropods.csv", header=T)

## The following code is useful for the newest version of tRophicPosition
### In this data set we have a combination of Study and Location variables, which need to be concatenated. Here we introduce mutate and arrange functions from the package dplyr and the pipe operator %>% (from magrittr but loaded through dplyr).

pterosTP <- pterosTP %>% mutate(Community = paste(Study,"-", Location, sep = ""))
pterosTP <- pterosTP %>% arrange(NS)

PteroList <- extractIsotopeData(pterosTP, b1 = "Pelagic_BL", b2 = "Pelagic_BL",
                                baselineColumn = "FG", consumersColumn = "Spp",
                                groupsColumn = "Community", 
                                d13C = "d13C", d15N = "d15N")

#use the str function to check that all went well during the use of the extractIsotopeData() function:
str(PteroList)

for(community in PteroList){
  print(summary(community))
  plot(community)
}

Ptero_models <- multiSpeciesTP(PteroList, model = "twoBaselinesFull",
                               n.adapt = 10000, n.iter = 10000,
                               burnin = 10000, n.chains = 5, print = FALSE)

# By default the mode is used in both trophic position and alpha plots
credibilityIntervals(Ptero_models$df, x = "group", xlab ="Community")

# Median instead of the mode,
# just add y1 and y2 as arguments
credibilityIntervals(Ptero_models$df, x = "group", xlab ="Community", 
                     y1 = "median", y2 = "alpha.median")

# To get a numerical summary
sapply(Ptero_models$"TPs", quantile, probs = c(0.025, 0.5, 0.975)) %>% round(3)

# To get the mode
getPosteriorMode(Ptero_models$"TPs")

# The following code was worked for an earlier version (0.7.0) of this package
## From here, obtain trophic position estimates for each species separately

## For Clione limacina 
clione <- loadIsotopeData(pterosTP,
                          species = "C. limacina",
                          b1 = "Benthic_BL",
                          b2 = "Pelagic_BL",
                          community = "oceanic",
                          speciesColumn = "FG",
                          baselineColumn = "FG",
                          communityColumn = "Location")
plot(clione, b1 = "Benthic_BL", b2 = "Pelagic_BL")

## The Bayesian model
model.string <- jagsBayesianModel(model = "twoBaselinesFull")
model <- TPmodel(data = clione, model.string = model.string,
                 n.adapt = 20000, n.chains = 2)
posterior.samples <- posteriorTP(model = model, n.iter = 20000,
                                 variable.names = c("TP", "muDeltaN"))
summary(posterior.samples)
posterior.combined <- coda::mcmc(do.call(rbind, posterior.samples))

## Calculate the mode 
getPosteriorMode(posterior.combined)
MCMCvis::MCMCtrace(posterior.samples)

# MCMCplot visualize posterior distributions using caterpillar plots
MCMCvis::MCMCplot(posterior.samples, ylim = c(0,5), horiz = FALSE)
mt <- as.matrix(posterior.samples)
clioneTP <- mt[,1]
head(clioneTP)

## For Spongiobranchaea australis
spongio <- loadIsotopeData(pterosTP,
                           species = "S. australis",
                           b1 = "Benthic_BL",
                           b2 = "Pelagic_BL",
                           community = "oceanic",
                           speciesColumn = "FG",
                           baselineColumn = "FG",
                           communityColumn = "Location")
plot(spongio, b1 = "Benthic_BL", b2 = "Pelagic_BL")

## The Bayesian model
model.string2 <- jagsBayesianModel(model = "twoBaselinesFull")
model2 <- TPmodel(data = spongio, model.string = model.string2,
                  n.adapt = 20000, n.chains = 2)
posterior.samples2 <- posteriorTP(model = model2, n.iter = 20000,
                                  variable.names = c("TP", "muDeltaN"))
summary(posterior.samples2)
posterior.combined2 <- coda::mcmc(do.call(rbind, posterior.samples2))

## Calculate the mode
getPosteriorMode(posterior.combined2)
MCMCvis::MCMCtrace(posterior.samples2)

# MCMCplot visualize posterior distributions using caterpillar plots
MCMCvis::MCMCplot(posterior.samples2, ylim = c(0,5), horiz = FALSE)
mt2 <- as.matrix(posterior.samples2)
spongioTP <- mt2[,1]
head(spongioTP)

## For Clio pyramidata
clioTP <- read.csv("trophicposition pteropods.csv", header=T)
clio <- loadIsotopeData(pterosTP,
                        species = "C. pyramidata",
                        b1 = "Benthic_BL",
                        b2 = "Pelagic_BL",
                        community = "oceanic",
                        speciesColumn = "FG",
                        baselineColumn = "FG",
                        communityColumn = "Location")
plot(clio, b1 = "Benthic_BL", b2 = "Pelagic_BL")

## The Bayesian model
model.string3 <- jagsBayesianModel(model = "twoBaselinesFull")
model3 <- TPmodel(data = clio, model.string = model.string3,
                  n.adapt = 20000, n.chains = 2)
posterior.samples3 <- posteriorTP(model = model3, n.iter = 20000,
                                  variable.names = c("TP", "muDeltaN"))
summary(posterior.samples3)
posterior.combined3 <- coda::mcmc(do.call(rbind, posterior.samples3))

## Calculate the mode
getPosteriorMode(posterior.combined3)
MCMCvis::MCMCtrace(posterior.samples3)

## MCMCplot visualize posterior distributions using caterpillar plots
MCMCvis::MCMCplot(posterior.samples3, ylim = c(0,5), horiz = FALSE)
mt3 <- as.matrix(posterior.samples3)
clioTP <- mt3[,1]
head(clioTP)

## Combine all into one matrix
pteropodsTP <- cbind(clioTP, clioneTP, spongioTP)
head(pteropodsTP)

## Melt
TPall <- melt(pteropodsTP, id.var = c('TP'), variable.name = 'species')
head(TPall)

## Change column headers
colnames(TPall) <- c("samples", "species", "TP")
head(TPall)

### I now have three objects containing posterior sampled TP values for the three spp: clioTP, clioneTP, and spongioTP; I want to combine them into a joy plot (aka ridge plot)

## Read image files
icon1 <- readPNG("Clio pyramidata_yellow.png")
icon2 <- readPNG("Clione limacina antarctica_grey.png")
icon3 <- readPNG("Spongiobranchaea australis_purple.png")

## Add to plot
a1 <- annotation_raster(icon1, ymin = 1, ymax = 1.6, xmin = 3.3, xmax=3.5) #clio
a2 <- annotation_raster(icon2, ymin = 2, ymax = 2.6, xmin = 3.6, xmax = 3.75) #clione
a3 <- annotation_raster(icon3, ymin = 3, ymax = 3.6, xmin = 2.85, xmax = 3) #spongio

## Create colour palette
clrs<-c(clioneTP = "#0E0D0D",clioTP = "#CCBA71",spongioTP = "#9985A5")

## Create ridge plot
TPjoy <- ggplot(TPall, aes(x = TP, y = species,
                           fill = species)) + geom_density_ridges(alpha = 0.7)

TPjoy1 <- TPjoy + a1 + a2 + a3
TPjoy1

## Add images to plot
TPjoy1 + scale_y_discrete("species", 
                          labels = expression("spongioTP"=italic("S. australis"), 
                                              "clioTP" = italic("C. pyramidata"), 
                                              "clioneTP" = italic("C. limacina"))) +
  labs(x = "trophic position") +
  scale_fill_cyclical(values = c("#CCBA71","#0E0D0D","#9985A5")) +
  xlim(2.5,4.5) +
  theme_bw() + theme(panel.grid = element_blank(), 
                     panel.border = element_blank(),
                     axis.title.y = element_blank(),
                     panel.grid.major.y = element_line(size = 0.5, 
                                                       color = "black")) 

## Figures A1 to A4
## The following script provides code for regression analyses testing d13C and d15N against spatial and temporal variables

library(ggplot2)
library(gridExtra)

## Read in data
SIAplots <- read.csv("dataframe for corrected isotopes plots.csv", header = TRUE)

## Omit NAs
SIAplots <- na.omit(SIAplots)

## Least squares regressions
my.model <- lm(d15N~d13C, data = SIAplots)
my.model 
summary(my.model)

## Rename each species for subsetting
clio <- subset(SIAplots, species=="C. pyramidata")
clione <- subset(SIAplots, species=="C. antarctica")
spongio <- subset(SIAplots, species=="S. australis")

## Want to combine day, month, and year columns to create one date column
SIAdate <- as.Date(with(SIAplots, paste(year, month, day, sep = "-")), "%Y-%m-%d")
SIAdate

SIAplots$Date <- paste( month.abb[SIAplots$month], SIAplots$day, sep = " " )

## Create some functions
lat <- SIAplots$lat
long <- SIAplots$lon

## Function to calculate mean and sd for each group
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

## Customise colours
cols <- c("C. pyramidata" = "#CCBA71",
          "C. antarctica" = "#0E0D0D",
          "S. australis" = "#9985A5")

## Customise labels
species_names <- c('C. pyramidata' = "Clio pyramidata",
                   'C. antarctica' = "Clione limacina antarctica",
                   'S. australis' = "Spongiobranchaea australis")

## ANOVAs
## All species ANOVA
anovaC <- aov(d13C~species + month + lat2 + long + depth2, data=SIAplots)
summary(anovaC) #species and month

anovaN <-aov(d15N~species + month + lat2 + long + depth2, data=SIAplots)
summary(anovaN) #species, latitude and longitude

## Per species
### Clio pyramidata
clioanovaC <- aov(d13C~month + lat2 + long + depth2, data=clio)
summary(clioanovaC) #month and longitude

clioanovaN <- aov(d15N~month + lat2 + long + depth2, data=clio)
summary(clioanovaN) #latitude

### Clione limacina antarctica
clioneanovaC <- aov(d13C~month + lat2 + long + depth2, data=clione)
summary(clioneanovaC) #month

clioneanovaN <- aov(d15N~month + lat2 + long + depth2, data=clione)
summary(clioneanovaN) #residuals

### Spongiobranchaea australis
spongioanovaC <- aov(d13C~lat2 + long + depth2, data=spongio)
summary(spongioanovaC) #month

spongioanovaN <- aov(d15N~lat2 + long + depth2, data=spongio)
summary(spongioanovaN) 

## Linear regression (isotopes ~ date) grouped by species 

### Clio pyramidata
clio_13C_date <- lm(d13C ~ day*month, data = clio)
summary(clio_13C_date) #R^2 = 0.21, p < 0.05
clio_15N_date <- lm(d15N ~ day*month, data = clio)
summary(clio_15N_date) #R^2 = 0.06, p < 0.05

### Clione limacina antarctica
clione_13C_date <- lm(d13C ~ day*month, data = clione)
summary(clione_13C_date) #R^2 = 0.50, p < 0.05
clione_15N_date <- lm(d15N ~ day*month, data = clione)
summary(clione_15N_date) #R^2 = 0.34, p < 0.05

### Spongiobranchaea australis
spongio_13C_date <- lm(d13C ~ day*month, data = spongio)
summary(spongio_13C_date) #R^2 = 0.11, p = 0.35
spongio_15N_date <- lm(d15N ~ day*month, data = spongio)
summary(spongio_15N_date) #R^2 = -0.39, p = 0.85

## Plot d13C vs date, grouped by species
d13Ctime2 <- ggplot(SIAplots, aes(x = SIAdate, y = d13C, color = species)) +
  geom_point() + 
  stat_smooth(method = "lm", formula = y ~ x, size = 1, se = FALSE) +
  labs(x = "",
       y = expression(paste(delta^13, "C (\u2030)", sep = ""))) +
  scale_colour_manual("species",
                      values = cols,
                      breaks = c("C. pyramidata",
                                 "C. antarctica",
                                 "S. australis"),
                      guide = FALSE) +guides(fill = FALSE) +
  facet_wrap(~species, labeller = as_labeller(species_names)) +
  theme_bw() +
  theme(strip.text = element_text(face = "italic")) 
d13Ctime2

## Plot d15N vs date, grouped by species 
d15Ntime2 <- ggplot(SIAplots, aes(x = SIAdate, y = d15N, color = species)) +
  geom_point() + 
  stat_smooth(method = "lm", formula = y ~ x, size = 1, se = FALSE) +
  labs(x = "Sampling dates",
       y = expression(paste(delta^15, "N (\u2030)", sep = ""))) +
  scale_colour_manual("species",
                      values = cols,
                      breaks = c("C. pyramidata",
                                 "C. antarctica",
                                 "S. australis"),
                      guide = FALSE) +guides(fill = FALSE) +
  facet_wrap(~species, labeller = as_labeller(species_names)) +
  theme_bw() +
  theme(strip.text = element_text(face = "italic"))
d15Ntime2

## Combine plots
grid.arrange(d13Ctime2, d15Ntime2, nrow = 2)

## Linear regression (isotopes ~ latitude) grouped by species 

### Clio pyramidata
clio_13C_lat <- lm(d13C ~ lat, data = clio)
summary(clio_13C_lat) #R^2 = -0.0001, p = 0.32
clio_15N_lat <- lm(d15N ~ lat, data = clio)
summary(clio_15N_lat) #R^2 = 0.04, p < 0.05

### Clione limacina antarctica
clione_13C_lat <- lm(d13C ~ lat, data = clione)
summary(clione_13C_lat) #R^2 = 0.03, p = 0.24
clione_15N_lat <- lm(d15N ~ lat, data = clione)
summary(clione_15N_lat) #R^2 = 0.02, p = 0.27

### Spongiobranchaea australis
spongio_13C_lat <- lm(d13C ~ lat, data = spongio)
summary(spongio_13C_lat) #R^2 = -0.14, p = 0.62
spongio_15N_lat <- lm(d15N ~ lat, data = spongio)
summary(spongio_15N_lat) #R^2 = -0.19, p = 0.84

## Plot d13C vs latitude, grouped by species
d13Clat2 <- ggplot(SIAplots, aes(x = lat, y = d13C, color = species)) + 
  geom_point() +
  stat_smooth(method = "lm", formula = y ~ x, size = 1, se = FALSE) +
  labs(x = "",
       y = expression(paste(delta^13, "C (\u2030)", sep = ""))) +
  scale_colour_manual("species",
                      values = cols,
                      breaks = c("C. pyramidata",
                                 "C. antarctica",
                                 "S. australis"),
                      guide = FALSE) +
  guides(fill = FALSE) +
  facet_wrap(~species, labeller = as_labeller(species_names)) +
  theme_bw() +
  theme(strip.text = element_text(face = "italic"))
d13Clat2

## Plot d13C vs latitude, grouped by species
d15Nlat2 <- ggplot(SIAplots, aes(x = lat, y = d15N, color = species)) + 
  geom_point() +
  stat_smooth(method = "lm", formula = y ~ x, size = 1, se = FALSE) +
  labs(x = "Latitude",
       y = expression(paste(delta^15, "N (\u2030)", sep = ""))) +
  scale_colour_manual("species",
                      values = cols,
                      breaks = c("C. pyramidata",
                                 "C. antarctica",
                                 "S. australis"),
                      guide = FALSE) +
  guides(fill = FALSE) +
  facet_wrap(~species, labeller = as_labeller(species_names)) +
  theme_bw() +
  theme(strip.text = element_text(face = "italic"))
d15Nlat2

## Combine plots
grid.arrange(d13Clat2, d15Nlat2, nrow = 2)



## Linear regression (isotopes ~ longitude) grouped by species 

### Clio pyramidata
clio_13C_lon <- lm(d13C ~ long, data = clio)
summary(clio_13C_lon) #R^2 = 0.06, p < 0.05
clio_15N_lon <- lm(d15N ~ long, data = clio)
summary(clio_15N_lon) #R^2 = 0.005, p = 0.20

### Clione limacina antarctica
clione_13C_lon <- lm(d13C ~ long, data = clione)
summary(clione_13C_lon) #R^2 = 0.20, p < 0.05
clione_15N_lon <- lm(d15N ~ long, data = clione)
summary(clione_15N_lon) #R^2 = -0.001, p = 0.34

### Spongiobranchaea australis
spongio_13C_lon <- lm(d13C ~ long, data = spongio)
summary(spongio_13C_lon) #R^2 = 0.07, p = 0.28
spongio_15N_lon <- lm(d15N ~ long, data = spongio)
summary(spongio_15N_lon) #R^2 = -0.11, p = 0.55

## Plot d13C vs longitude, grouped by species
d13Clong2 <- ggplot(SIAplots, aes(x = long, y = d13C, color = species)) + 
  geom_point() +
  stat_smooth(method = "lm", formula = y ~ x, size = 1, se = FALSE) +
  labs(x = " ",
       y = expression(paste(delta^13, "C (\u2030)", sep = ""))) +
  scale_colour_manual("species",
                      values = cols,
                      breaks = c("C. pyramidata",
                                 "C. antarctica",
                                 "S. australis"),
                      guide = FALSE) +
  guides(fill = FALSE) +
  facet_wrap(~species, labeller = as_labeller(species_names)) +
  theme_bw() +
  theme(strip.text = element_text(face = "italic"))
d13Clong2

## Plot d15N vs longitude, grouped by species
d15Nlong2 <- ggplot(SIAplots, aes(x = long, y = d15N, color = species)) + 
  geom_point() +
  stat_smooth(method = "lm", formula = y ~ x, size = 1, se = FALSE) +
  labs(x = "Longitude",
       y = expression(paste(delta^15, "N (\u2030)", sep = ""))) +
  scale_colour_manual("species",
                      values = cols,
                      breaks = c("C. pyramidata",
                                 "C. antarctica",
                                 "S. australis"),
                      guide = FALSE) +
  guides(fill = FALSE) +
  facet_wrap(~species, labeller = as_labeller(species_names)) +
  theme_bw() +
  theme(strip.text = element_text(face = "italic"))
d15Nlong2

## Combine plots
grid.arrange(d13Clong2, d15Nlong2, nrow = 2)


#Linear regression (isotopes ~ depth) grouped by species 

### Clio pyramidata
clio_13C_depth <- lm(d13C ~ depth, data = clio)
summary(clio_13C_depth) #R^2 = 0.006, p = 0.17
clio_15N_depth <- lm(d15N ~ depth, data = clio)
summary(clio_15N_depth) #R^2 = -0.006, p = 0.67

### Clione limacina antarctica
clione_13C_depth <- lm(d13C ~ depth, data = clione)
summary(clione_13C_depth) #R^2 = -0.05, p = 0.62
clione_15N_depth <- lm(d15N ~ depth, data = clione)
summary(clione_15N_depth) #R^2 = 0.14, p = 0.07

### Spongiobranchaea australis
spongio_13C_depth <- lm(d13C ~ depth, data = spongio)
summary(spongio_13C_depth) #R^2 = -0.13, p = 0.60
spongio_15N_depth <- lm(d15N ~ depth, data = spongio)
summary(spongio_15N_depth) #R^2 = -0.13, p = 0.59

## Plot d13C vs depth, grouped by species
d13Cdepth2 <- ggplot(SIAplots, aes(x = depth, y = d13C, color = species)) + 
  geom_point() +
  stat_smooth(method = "lm", formula = y ~ x, size = 1, se = FALSE) +
  labs(x = " ",
       y = expression(paste(delta^13, "C (\u2030)", sep = ""))) +
  scale_colour_manual("species",
                      values = cols,
                      breaks = c("C. pyramidata",
                                 "C. antarctica",
                                 "S. australis"),
                      guide = FALSE) +
  guides(fill = FALSE) +
  facet_wrap(~species, labeller = as_labeller(species_names)) +
  theme_bw() +
  theme(strip.text = element_text(face = "italic"))
d13Cdepth2

## Plot d15N vs depth, grouped by species
d15Ndepth2 <- ggplot(SIAplots, aes(x = depth, y = d15N, color = species)) + 
  geom_point() +
  stat_smooth(method = "lm", formula = y ~ x, size = 1, se = FALSE) +
  labs(x = "Depth (m)",
       y = expression(paste(delta^15, "N (\u2030)", sep = ""))) +
  scale_colour_manual("species",
                      values = cols,
                      breaks = c("C. pyramidata",
                                 "C. antarctica",
                                 "S. australis"),
                      guide = FALSE) +
  guides(fill = FALSE) +
  facet_wrap(~species, labeller = as_labeller(species_names)) +
  theme_bw() +
  theme(strip.text = element_text(face = "italic"))
d15Ndepth2

## Combine plots
grid.arrange(d13Cdepth2, d15Ndepth2, nrow = 2)
