##################################################################################################################
## Project: Prepare available wild fire perimeter data for use in studying potential outcomes
##
## Script purpose: Read weather variables created by FireFamilyplus v5 software to feed climate conditioning
##                 process
##
## Date: 25th March 2020
## Author: JR
##################################################################################################################

# Required testing many libraries. Be selective in initializing libraries
library(rlang) 
library(maptools) 
library(GISTools) 
library(plyr)
library(dplyr)
library(lwgeom)
library(sf)
library(sp)
library(stringr)
library(smoothr)
library(geosphere)
library(rasterVis) 
library(gdalUtils)
library(rgdal)
library(raster) 
library(rgeos)
library(velox)
library(parallel)
library(data.table)
library(geojsonio)
library(cleangeo)

##################################################################################################################
## Section: Set global parameters
##################################################################################################################

# Projection string for meter projection US wide
AEAProj <- '+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-110
+x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m'

# Sheffs preferred output projections
laea_proj4 <- '+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs'

# Test for valid GDAL install
gdal_setInstallation()
valid_install <- !is.null(getOption('gdalUtils_gdalPath'))
if(require(raster) && require(rgdal) && valid_install)
{
  print('TRUE')
}
#gdal_setInstallation(ignore.full_scan=FALSE)

# Set raster options
rasterOptions(maxmemory = 5.0e+09)

#if("package:adehabitatHR" %in% search()) detach("package:adehabitatHR", unload=TRUE)

memory.limit()
# set max memory usage 
memory.size(max=5.0e+09)

##################################################################################################################
## Section: Set global directory locations
##################################################################################################################

# Check directory
getwd()
wd <- getwd()
print(paste0('Current working dir: ', wd))

# Temporary raster location
raster::rasterOptions(tmpdir = wd)

# Get the data directory from readtext
DATA_DIR <- ('Data') 
OUTPUT_DIR <- ('Output') 

##################################################################################################################
## Section: Load starting wild fire perimeter data
##
## Data consiSts of three data sets in shapefile format:
##
## 1. Perimeters 1984 thru 2017 inclusive, from the MTBS. (https://www.mtbs.gov/)
##    Data called, 'mtbs_perims_DD.xxx' Data represents final event perimeters.
##
## 2. 2018 and 2019 perimeter data from GEOMACS. (https://www.geomac.gov/)
##    Data called, '2018_perimeters_dd83.xxx' and '2019_perimeters_dd83.xxx'
##    Date represents temporal spread of events. Multiple perimeters per event.
##    Key to using the data is to create "Final" event perimeters.
##
## The aim is to load and combine data sets ellimnating data not covered by the current risQ geographic domain.
## As such none coterminous US data will be removed. The attribute data contained by each data set will be unified
## likely reducing the attribute fields in the final data set
##
## All data was downlaoded 25th March 2020
##
##################################################################################################################

# Set input directory
initial_path <- file.path(paste0('./',DATA_DIR))

###################################### MTBS DATA LOAD ###################################### 

# Read MTBS shapefile
mtbs_shp <- readOGR(paste0(initial_path,'/','mtbs_perims_DD.shp')) 
mtbs_shp <- spTransform(mtbs_shp, CRS('+init=epsg:4326'))

# Inspect attribute data (22967 perimeters with 7 attributes, 6 factor & 1 numeric)
# Factor fields require conversion to character in order to work with or join on them
str(mtbs_shp@data)
head(mtbs_shp@data)
levels(mtbs_shp@data$Fire_Type)
#plot(mtbs_shp)

# Data set contains AK, HI and PR from viewing the plot - require to be removed, state code is embedded in 'Fire_ID'
# Creat efield containing state abbreviations
mtbs_shp@data$NEWFIREID <-str_trim(as.character(mtbs_shp@data$Fire_ID))
mtbs_shp@data$STATE <- substring((as.character(mtbs_shp@data$NEWFIREID)), 1, 2)  
head(mtbs_shp@data)
tail(mtbs_shp@data)

# Get unique list of abbreviated state names
uniq.sts.df <- unique(unlist(sapply(mtbs_shp@data$STATE, function(x) unique(x[]))))

# R has built in state.abb and state.nme function to create list of coterminous US state abbreviation to differnce with
# If in doubt refer to https://www.infoplease.com/us/postal-information/state-abbreviations-and-state-postal-codes
non.us.sts <- setdiff(unlist(uniq.sts.df),state.abb)

# Require to remove just Alazka, Hawaii and Puerto Rico - AK, HI & PR
mtbs.us.ml.only <- subset(mtbs_shp, !(STATE %in% c('AK', 'HI', 'PR')))
#plot(mtbs.us.ml.only)

# Test for uniqueness of fire ids - Distinct event perimiters are unique uni.mtbs.evnts = nrow(mtbs.us.ml.only.DT)
mtbs.us.ml.only.DT <- as.data.table(mtbs.us.ml.only@data)
uni.mtbs.evnts <- unique(mtbs.us.ml.only.DT, by = c('Fire_ID')) 

# Save current output
writeOGR(obj=mtbs.us.ml.only, dsn=final_path, layer='mtbs_perims_DD_Clnd', driver="ESRI Shapefile", overwrite=TRUE)

###################################### MTBS DATA LOAD ###################################### 

###################################### GEOMACS DATA LOAD ###################################

# Read GEOMACS shapefiles
toMatch <- c('dd83','2019','2018')
shp.list <- list.files(path=initial_path, pattern=glob2rx('*shp*'), full.names=T,recursive=F)
fls.nded <- unique(grep(paste(toMatch,collapse="|"), shp.list, value=TRUE))
fls.nded <- exclude(fls.nded,'xml')

# Attribute feature columns don't match, there require to read files individually
#do.call(rbind, lapply(fls.nded, rgdal::readOGR))
GEO18_shp <- readOGR(fls.nded[1]) 
GEO18_shp <- spTransform(GEO18_shp, CRS('+init=epsg:4269'))
GEO19_shp <- readOGR(fls.nded[2]) 
GEO19_shp <- spTransform(GEO19_shp, CRS('+init=epsg:4269'))
#plot(GEO18_shp)
#plot(GEO19_shp, add = TRUE)

# Examine difference in fireld names
difs <- setdiff(names(GEO18_shp@data),names(GEO19_shp@data))
difs2 <- setdiff(names(GEO19_shp@data),names(GEO18_shp@data))

# Rename fields
names(GEO18_shp@data)[names(GEO18_shp@data)=='temp'] <- 'objectid'
names(GEO19_shp@data)[names(GEO19_shp@data)=='st_area_sh'] <- 'shape_Area'
names(GEO19_shp@data)[names(GEO19_shp@data)=='st_length_'] <- 'shape_Leng'

# Append 2018 ith2019 data
GEOMgd <- rbind(GEO18_shp, GEO19_shp, makeUniqueIDs = TRUE)  
#plot(GEOMgd)

# Add row ID to attributes
GEOMgd@data$idin <- attr(GEOMgd@data, 'row.names') # Returns integers
head(GEOMgd@data)

# Get unique list of abbreviated state names
uniq.stsl.df <- unique(unlist(sapply(GEOMgd@data$state, function(x) unique(x[]))))

# None coterminous US states
non.us.stsl <- setdiff(unlist(uniq.stsl.df),state.abb)

# Require to remove just Alazka, Hawaii and Puerto Rico - AK, HI & PR
GEOMgd.only <- subset(GEOMgd, !(state %in% c('AK', 'HI')))
levels(droplevels(GEOMgd.only$state))

# Save current output
writeOGR(obj=GEOMgd.only, dsn=final_path, layer='wfirepmtrs2018-2019', driver="ESRI Shapefile", overwrite=TRUE)

# Test for uniqueness of fire ids - Distinct event perimiters are not unique uni.GEO.evnts < nrow(GEOMgd.only.DT)
GEOMgd.only.DT <- as.data.table(GEOMgd.only@data)

setorder(GEOMgd.only.DT, uniquefire, latest)
GEOMgd.only.DT[ , `:=`( COUNT = .N , IDX = 1:.N ) , by = c('uniquefire') ]

uni.GEO.evnts <- unique(GEOMgd.only.DT, by = c('uniquefire'))  
uni.y.ltst <- GEOMgd.only.DT[latest == 'Y']
uni.y.ltst.ts <- unique(uni.y.ltst, by = c('uniquefire')) 

# Appears selecting latest = Y, provides most unique result, 6 fires are dropped. 
GEOMgd.fnl01 <- GEOMgd.only[GEOMgd.only$latest == 'Y',]
levels(droplevels(GEOMgd.fnl01$latest))

# Reproject data to match projection of MTBS data
GEOMgd.fnl.rpj <- spTransform(GEOMgd.fnl01,crs(mtbs.us.ml.only))

# Save current output
writeOGR(obj=GEOMgd.fnl.rpj, dsn=final_path, layer='wfirepmtrs2018-2019_LatestPmters2', driver="ESRI Shapefile", overwrite=TRUE)

###################################### GEOMACS DATA LOAD ###################################

##################################################################################################################
## Section: Merge MTBS data with the GEOMACS combined data 
##################################################################################################################

# Compare field names from each data source - joining will require normalizing names and column count
names(mtbs.us.ml.only@data)
names(GEOMgd.fnl.rpj@data)

# Basic data required:
# 1. Unique fire id (Fire_ID)
# 2. Fire name      (Fire_Name)
# 3. Year           (Year)
# 4. Fire type      (Fire_Type)
# 5. Size (Acres)   (Acres)
# 6. State          (STATE)

# Cleen MTBS data names - require to drop columns StartMonth StartDay & NEWFIREID
dropcs <- c('StartMonth','StartDay','NEWFIREID') # list of col names
mtbs.fnl <- mtbs.us.ml.only[,!(names(mtbs.us.ml.only) %in% dropcs)]

# Cleen GEOMACS data names - laborious way
names(GEOMgd.fnl.rpj@data)[names(GEOMgd.fnl.rpj@data)=='uniquefire'] <- 'Fire_ID'
names(GEOMgd.fnl.rpj@data)[names(GEOMgd.fnl.rpj@data)=='incidentna'] <- 'Fire_Name'
names(GEOMgd.fnl.rpj@data)[names(GEOMgd.fnl.rpj@data)=='fireyear'] <- 'Year'
names(GEOMgd.fnl.rpj@data)[names(GEOMgd.fnl.rpj@data)=='gisacres'] <- 'Acres'
names(GEOMgd.fnl.rpj@data)[names(GEOMgd.fnl.rpj@data)=='state'] <- 'STATE'
names(GEOMgd.fnl.rpj@data)[names(GEOMgd.fnl.rpj@data)=='firecode'] <- 'Fire_Type'

# Keep GEOMACS attribute fields not required
keepcs <- c('Fire_ID','Fire_Name','Year','Acres','STATE','Fire_Type')
GEOMgd.fnl <- GEOMgd.fnl.rpj[,(names(GEOMgd.fnl.rpj) %in% keepcs)]

# Append MTBS and GEOMACS data
wfirepmtrs <- rbind(mtbs.fnl, GEOMgd.fnl, makeUniqueIDs = TRUE)  
#plot(wfirepmtrs)

# Set output directory
final_path <- file.path(paste0('./',OUTPUT_DIR))

# Save current output
writeOGR(obj=wfirepmtrs, dsn=final_path, layer='wfirepmtrs1984-2019', driver="ESRI Shapefile", overwrite=TRUE)

##################################################################################################################
## Section: Drop fires based on cause and size
##################################################################################################################

# Re-read wild fire perimeters file
wfirepmtrs_shp <- readOGR(paste0(final_path,'/','wfirepmtrs1984-2019.shp')) 
wfirepmtrs_shp <- spTransform(wfirepmtrs_shp, CRS('+init=epsg:4326'))

# Remove prescribed burns, likely most from MTBS data (Oklahoma bulls eye of fire activity, will be removed)
wfirepmtrs_shp <- subset(wfirepmtrs_shp, !(Fire_Type %like% c('Prescribed')))

# Reproject data to enable accurate polygon area calculation - output km2 converted acres 
wfirepmtrs.rpj <- spTransform(wfirepmtrs_shp,crs(laea_proj4))
wfirepmtrs.rpj@data$Acres <- raster::area(wfirepmtrs.rpj)*247.11/1000000.0
wfirepmtrs.DD <- spTransform(wfirepmtrs.rpj,CRS('+init=epsg:4326'))

# Filter events based on 300 Acre cutoff
wfirepmtrs.DD <- subset(wfirepmtrs.DD, Acres >= 50.0)

# Filter errors on Oklahoma - obvious prescribed burns coded as unknown
wfirepmtrs.DD <- subset(wfirepmtrs.DD, !(STATE %in% c('OK','KS') & Fire_Type %in% c('Unknown')))
#plot(wfirepmtrs.DD)

# Save current output
writeOGR(obj=wfirepmtrs.DD, dsn=final_path, layer='fltrdwfirepmtrs1984-2019', driver="ESRI Shapefile", overwrite=TRUE)

##################################################################################################################
## Section: Intersect final perimeters with WUI mask
##################################################################################################################

# Read WUI mask raster and select pixels value equal to 1
#WUI.msk <- raster::raster(paste0(initial_path,'/','WFUrbanMask.tif'))  
#WUI.msk[WUI.msk != 1] <- 0

# Reproject region extent - set type
#WUI.msk.ext <- projectExtent(WUI.msk, CRS('+init=epsg:4326'))
#res(WUI.msk.ext) <- 0.003
#WUI.msk.prj <- projectRaster(WUI.msk, WUI.msk.ext, method = 'ngb')
#plot(WUI.msk.prj)

# Set reprojected ratser to integer type
#WUI.msk.prj[]=as.integer(WUI.msk.prj[])
#storage.mode(WUI.msk.prj[])

# Save reprojected raster
#if (require(rgdal)) {
#  rf <- writeRaster(WUI.msk.prj, filename=paste0(final_path,'/','WFUrbanMask_DD.tif'), 
#                    format='GTiff', datatype='INT4S', overwrite=TRUE)
#}

# Re-read WUI mask raster 
WUI.msk <- raster::raster(paste0(final_path,'/','WFUrbanMask_DD.tif')) 

# Re-read, clean and write updated file as geojson
wfirepmtrs.DD <- readOGR(paste0(final_path,'/','fltrdwfirepmtrs1984-2019.shp')) 
wfirepmtrs.DD@data$idin <- attr(wfirepmtrs.DD@data, 'row.names') # Returns integers
#plot(wfirepmtrs.DD)

# Overlay wildfire perimeters with urban mask - most efficient way, count pixels per perimeter and drop 
# perimeters with few pixels
ras.vx <- velox(WUI.msk)

# Count WUI pixels within each fire perimeter
result <- mclapply(seq_along(1), function(x){
  q <- ras.vx$crop(wfirepmtrs.DD);ras.vx$extract(wfirepmtrs.DD, fun=function(t) sum(t,na.rm=T))
})
result.dt <- as.data.table(Reduce(rbind, result))
result.dt[ , `:=`( idin = 0:(.N-1) )]
setnames(result.dt, 'V1', 'COUNT')
polylst <- result.dt[COUNT > 0][, idin]

# Hand test polygons if necessary
test2 <- wfirepmtrs.DD[5,'idin']
bbin <- st_bbox(test2)
xminbb = bbin[[1]]
xmaxbb = bbin[[3]]
yminbb = bbin[[2]]
ymaxbb = bbin[[4]]
xlimin = c(xminbb, xmaxbb) 
ylimin = c(yminbb, ymaxbb) 
plot(WUI.msk, xlim = xlimin, ylim = ylimin)
plot(test2, add = TRUE)

# Subset final wildfire perimeters based on WUI pixel count greater than zero
wfirepmtrs.DD.fltdrd <- subset(wfirepmtrs.DD, (idin %in% polylst))
#plot(wfirepmtrs.DD.fltdrd)
#plot.new()

# Read attributes
wfirepmtrs.DD.fltdrd.df <- as.data.frame(wfirepmtrs.DD.fltdrd@data)

# Ensure poorly formed polygons are fixed
wfirepmtrs.DD.fltdrd.smpl <- gSimplify(wfirepmtrs.DD.fltdrd, tol = 0.00001, topologyPreserve=TRUE)

# Create a spatial polygon data frame 
wfirepmtrs.smpl <- SpatialPolygonsDataFrame(wfirepmtrs.DD.fltdrd.smpl, wfirepmtrs.DD.fltdrd.df)

# Ensure complete polygons
wfirepmtrs.DD.fltdrd <- spTransform(wfirepmtrs.smpl,crs(laea_proj4))
wfirepmtrs.DD.fltdrd <- gBuffer(wfirepmtrs.DD.fltdrd, byid=TRUE, width=0)
wfirepmtrs.DD.fltdrd <- spTransform(wfirepmtrs.DD.fltdrd,CRS('+init=epsg:4326'))

# Test for bad polygons
sum(gIsValid(wfirepmtrs.DD.fltdrd, byid=TRUE)==FALSE)

# Miscellaneous testing
#testout <- st_as_sfc(wfirepmtrs.DD.fltdrd)
#testout <- lwgeom::st_make_valid(testout)
#testout.vd <- st_is_valid(testout)
#which(testout.vd == FALSE, arr.ind = TRUE)
#sp.clean <- clgeo_Clean(wfirepmtrs.DD.fltdrd) 
#report.clean <- clgeo_CollectionReport(sp.clean) 
#clgeo_SummaryReport(report.clean) 

# Save current output
writeOGR(obj=wfirepmtrs.DD.fltdrd, dsn=final_path, layer='final_wfirepmtrs1984-2019', driver="ESRI Shapefile", overwrite=TRUE)

# Output geojson final files post validate and fix
geojson_write(wfirepmtrs.DD.fltdrd, lat = NULL, lon = NULL, geometry = 'polygon',
              group = NULL, file = paste0(final_path,'/','finalwfirepmtrs1984-2019_clnd.geojson'), overwrite = TRUE,
              precision = 6, convert_wgs84 = FALSE, crs = CRS('+init=epsg:4326'))

##################################################################################################################
## Section: Read accumulated exposures data aggregated through BQuery for each polygon
##################################################################################################################

# Read exposures data
final_path <- file.path(paste0('./',DATA_DIR))
expofile <- paste0(final_path,'/','wildfire-perimeter-totals.csv')
expo.DT <- read.csv(expofile)
setDT(expo.DT)

# Drop obsolete columns
expo.DT[ ,`:=`(boi_type = NULL ,radius = NULL, mode = NULL)] 

# Re-read wild fire perimeters file exposures aggregations were based on
final_path <- file.path(paste0('./',OUTPUT_DIR))
wfirepmtrs_shp <- readOGR(paste0(final_path,'/','final_wfirepmtrs1984-2019.shp')) 
wfirepmtrs_shp <- spTransform(wfirepmtrs_shp, CRS('+init=epsg:4326'))

# Get attribute data from boundary file
wfirepmtrs_shp.DT <- wfirepmtrs_shp@data
setDT(wfirepmtrs_shp.DT)

# Merge exposures data with wildfire perimeter attributes
setkey(wfirepmtrs_shp.DT, Fire_ID)
setkey(expo.DT, boi_source_link_id)
expo.DT <- merge(expo.DT, wfirepmtrs_shp.DT, by.x='boi_source_link_id', by.y='Fire_ID', all = FALSE)

# Drop obsolete columns
expo.DT[ ,`:=`(Fire_Type = NULL)] 

# Test data Extract Camp fire
expo.DT[Fire_Name == 'PARADISE']

##################################################################################################################
## Section: Define miscellaneous functions
##################################################################################################################

#Returns all items in a list that are not contained in toMatch
#toMatch can be a single item or a list of items
exclude <- function (theList, toMatch){
  return(setdiff(theList,include(theList,toMatch)))
}

#Returns all items in a list that ARE contained in toMatch
#toMatch can be a single item or a list of items
include <- function (theList, toMatch){
  matches <- unique (grep(paste(toMatch,collapse="|"), 
                          theList, value=TRUE))
  return(matches)
}

# Base R function to remove leading or trailing white spaces
trim <- function( x ) {
  gsub("(^[[:space:]]+|[[:space:]]+$)", "", x)
}

# Test if year is leap or not
is.leapyear <- function(year){
  return(((year %% 4 == 0) & (year %% 100 != 0)) | (year %% 400 == 0))
}

# Test sequence of numbers for consecutive differences of 1, can be changed
seqle <- function(x,incr=1) { 
  if(!is.numeric(x)) x <- as.numeric(x) 
  n <- length(x)  
  #y <- x[-1L] != x[-n] + incr 
  y <- abs(x[-1L] - x[-n] - incr) > .Machine$double.eps ^ 0.5
  i <- c(which(y|is.na(y)),n) 
  list(lengths = diff(c(0L,i)),
       values = x[head(c(0L,i)+1L,-1L)]) 
} 

##################################################################################################################
## Section: Book keeping - Clean memory close file connections
##################################################################################################################

# CLEAN MEMORY
rm(list = ls(all.names = TRUE))
raster::removeTmpFiles(h = 0)
flush.console()

##################################################################################################################
##                                                                                                              ##
##            Program end section Prepare wild fire perimeters for studying potential outcomes                  ##
##                                                                                                              ##
##################################################################################################################








