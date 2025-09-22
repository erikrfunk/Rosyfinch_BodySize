library("raster")
library("dismo")
library("rgeos")
library("rJava")

# If layers are in a single tif, run this first to convert to asc
r<-stack('~/PPM/sdm/bioclims/.tif')
r=raster::stack(list.files(paste0("~/finches/rosyfinches/BodySize_proj/climate_data/wc2.1_2.5m_bio/"), pattern = ".tif$", 
             full.names = T)) # Or run this if layers are multiple files
plot(r)
nlayers(r)
for(i in 1:nlayers(r)){
  band<-r[[i]]
  writeRaster(band,paste0('climate_data/',names(band),'.asc'),overwrite=TRUE)
}

# Then load the stack of asc rasters
clim_list <- list.files(paste0("~/finches/rosyfinches/BodySize_proj/climate_data/"), pattern = ".asc$", 
                        full.names = T) # Set for focal species only
clim <- raster::stack(clim_list)

occ = read.table("bcrf_winter_unique_coords.txt",header=T,sep="\t")
bounds = c(min(occ$Longitude)-2, -80,min(occ$Latitude)-2, max(occ$Latitude)+2) # Min and max observed lat/longs
e <- as(extent(bounds), 'SpatialPolygons') # bounding box for analyses 
crs(e) <- "+proj=longlat +datum=WGS84 +no_defs"
clim <- crop(clim, e)
#occ = occ[occ$lon > bounds[1] & occ$lon < bounds[2],]
#occ = occ[occ$lat > bounds[3] & occ$lat < bounds[4],]
dups = duplicated(occ[c("lat","lon")])
occ_unique = occ[!dups,]

coordinates(occ_unique) <- ~longitude + ltitude
cells <- cellFromXY(clim$wc2.1_2.5m_bio_6, occ_unique)
dups <- duplicated(cells)
occ_final <- occ_unique[!dups, ]
occ_buff <- buffer(occ_final, 4000)
studyArea <- crop(clim,extent(occ_buff))  

writeRaster(studyArea,
            # a series of names for output files
            filename=paste0("2040-2060_cropped/",names(studyArea),".asc"), 
            format="ascii", ## the output format
            bylayer=TRUE, ## this will save a series of layers
            overwrite=T)

set.seed(1) # For reproducibility
bg <- sampleRandom(x=clim,
                   size=min(10000,studyArea@nrows*studyArea@ncols),
                   na.rm=T, #removes the 'Not Applicable' points  
                   sp=T) # return spatial points
plot(studyArea[[2]])
plot(bg,add=T) 
plot(occ_final,add=T,col="black",pch = 1)

selected <- sample(1:length(occ_final), length(occ_final) * 0.75)
occ_train <- occ_final[selected, ]  # this is the selection to be used for model training
occ_test <- occ_final[-selected, ]  # this is the opposite of the selection which will be used for model testing
p <- raster::extract(clim, occ_train) # extracting env conditions as tabular data
p_test <- raster::extract(clim, occ_test)
a <- raster::extract(clim, bg)
pa <- c(rep(1, nrow(p)), rep(0, nrow(a))) # presence / absence
pder <- as.data.frame(rbind(p, a)) # climate predictors

mod <- maxent(x=pder, ## env conditions
              p=pa,   ## 1:presence or 0:absence
              path=paste0("/home/erik/PPM/sdm/tmpouts/"),
              args=c("responsecurves")) ## parameter specification

pred2010temp <- predict(mod, clim2010) #project to study area [raster]
plot(pred2010temp)

## Model evaluation
mod_eval_train <- dismo::evaluate(p = p, a = a, model = mod)
print(mod_eval_train)
mod_eval_test <- dismo::evaluate(p = p_test, a = a, model = mod)
print(mod_eval_test)

# Niche equivalency
occ_1 = read.table("DipMer_occ_final")
coordinates(occ_1) = ~lon+lat
occ_2 = read.table("DipAgi_occ_final")
coordinates(occ_2) = ~lon+lat
NE = nicheEquivalency(occs[[2]][,c("lon","lat")],occs[[1]][,c("lon","lat")],clim)

NO = nicheOverlap(pred1,pred2,stat="I")

