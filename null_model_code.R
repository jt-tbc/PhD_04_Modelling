
setwd("C:/working/Shapefiles/10m/Tables")
df_10m_1 <- read.csv("01_South_10mExtract_10m.csv")
df_10m_2 <- read.csv("02_North_10mExtract_10m.csv")
df_10m <- rbind(df_10m_1, df_10m_2)

setwd("C:/working/Shapefiles/25m/Tables")
df_25m <- read.csv("01_GT_10mExtract_25m.csv")

setwd("C:/working/Shapefiles/100m/Tables")
df_100m <- read.csv("01_GT_10mExtract_100m.csv")

# coral prevalence
coral.p <- length(which(df_25m$CORAL == 1)) / length(df_25m$CORAL)
coral.p2 <- length(which(df_100m$CORAL == 1)) / length(df_100m$CORAL)
coral.p3 <- length(which(df_10m$CORAL == 1)) / length(df_10m$CORAL)

# MACRO prevalence
MACRO.p <- length(which(df_25m$MACRO == 1)) / length(df_25m$MACRO)
MACRO.p2 <- length(which(df_100m$MACRO == 1)) / length(df_100m$MACRO)
MACRO.p3 <- length(which(df_10m$MACRO == 1)) / length(df_10m$MACRO)

# SPONGE prevalence
SPONGE.p <- length(which(df_25m$SPONGE == 1)) / length(df_25m$SPONGE)
SPONGE.p2 <- length(which(df_100m$SPONGE == 1)) / length(df_100m$SPONGE)
SPONGE.p3 <- length(which(df_10m$SPONGE == 1)) / length(df_10m$SPONGE)

# Load rasters


# Length of rasters
l.10m <- length()
l.50m <- length()
l.100m <- length()
l.250m <- length()

# Create null layers
null.10m_CORAL <- rep(NA,l.10m)
null.50m_CORAL <- rep(NA,l.50m)
null.100m_CORAL <- rep(NA,l.100m)
null.250m_CORAL <- rep(NA,l.250m)

null.10m_MACRO <- rep(NA,l.10m)
null.50m_MACRO <- rep(NA,l.50m)
null.100m_MACRO <- rep(NA,l.100m)
null.250m_MACRO <- rep(NA,l.250m)

null.10m_SPONGE <- rep(NA,l.10m)
null.50m_SPONGE <- rep(NA,l.50m)
null.100m_SPONGE <- rep(NA,l.100m)
null.250m_SPONGE <- rep(NA,l.250m)

for (i in 1:l.10m) null.10m_CORAL[i] <- rbinom(c(0,1),1,prob=c(1-coral.p,coral.p))
for (i in 1:l.50m) null.50m_CORAL[i] <- rbinom(c(0,1),1,prob=c(1-coral.p,coral.p))
for (i in 1:l.100m) null.100m_CORAL[i] <- rbinom(c(0,1),1,prob=c(1-coral.p,coral.p))
for (i in 1:l.250m) null.250m_CORAL[i] <- rbinom(c(0,1),1,prob=c(1-coral.p,coral.p))

for (i in 1:l.10m) null.10m_MACRO[i] <- rbinom(c(0,1),1,prob=c(1-MACRO.p,MACRO.p))
for (i in 1:l.50m) null.50m_MACRO[i] <- rbinom(c(0,1),1,prob=c(1-MACRO.p,MACRO.p))
for (i in 1:l.100m) null.100m_MACRO[i] <- rbinom(c(0,1),1,prob=c(1-MACRO.p,MACRO.p))
for (i in 1:l.250m) null.250m_MACRO[i] <- rbinom(c(0,1),1,prob=c(1-MACRO.p,MACRO.p))

for (i in 1:l.10m) null.10m_SPONGE[i] <- rbinom(c(0,1),1,prob=c(1-SPONGE.p,SPONGE.p))
for (i in 1:l.50m) null.50m_SPONGE[i] <- rbinom(c(0,1),1,prob=c(1-SPONGE.p,SPONGE.p))
for (i in 1:l.100m) null.100m_SPONGE[i] <- rbinom(c(0,1),1,prob=c(1-SPONGE.p,SPONGE.p))
for (i in 1:l.250m) null.250m_SPONGE[i] <- rbinom(c(0,1),1,prob=c(1-SPONGE.p,SPONGE.p))


# Compare to layers



# length of the raster
.d <- c(1,1,1,1,1,0,0,0,0,0,0,0)
# get .p from the data to determine probability of prescence, this is the prevelance
# use this to build a random layer
# you can then compare this surface to the generated layers
# wilcoxin test, non-parametric test, or paired t-test, check which test is appropriate
.p <- length(which(.d == 1)) / length(.d)
.pr <- rep(NA,12)
for (i in 1:12) .pr[i] <- rbinom(c(0,1),1,prob=c(1-.p,.p))
# generate p value about whether significantly different
