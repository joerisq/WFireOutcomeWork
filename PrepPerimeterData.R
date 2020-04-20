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

# Drop fires from ecah year not in current year
GEO18_shp.Tru18 <- GEO18_shp[GEO18_shp$fireyear == 2018,]
GEO19_shp.Tru19 <- GEO19_shp[GEO19_shp$fireyear == 2019,]

# Examine difference in fireld names
difs <- setdiff(names(GEO18_shp.Tru18@data),names(GEO19_shp.Tru19@data))
difs2 <- setdiff(names(GEO19_shp.Tru19@data),names(GEO18_shp.Tru18@data))

# Rename fields
names(GEO18_shp.Tru18@data)[names(GEO18_shp.Tru18@data)=='temp'] <- 'objectid'
names(GEO19_shp.Tru19@data)[names(GEO19_shp.Tru19@data)=='st_area_sh'] <- 'shape_Area'
names(GEO19_shp.Tru19@data)[names(GEO19_shp.Tru19@data)=='st_length_'] <- 'shape_Leng'

# Append 2018 with 2019 data
GEOMgd <- rbind(GEO18_shp.Tru18, GEO19_shp.Tru19, makeUniqueIDs = TRUE)  
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
writeOGR(obj=GEOMgd.fnl.rpj, dsn=final_path, layer='wfirepmtrs2018-2019_LatestPmters', driver="ESRI Shapefile", overwrite=TRUE)

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
levels(droplevels(wfirepmtrs_shp$Fire_Type))

# Reproject data to enable accurate polygon area calculation - output km2 converted to acres 
wfirepmtrs.rpj <- spTransform(wfirepmtrs_shp,crs(laea_proj4))
wfirepmtrs.rpj@data$Acres <- raster::area(wfirepmtrs.rpj)*247.11/1000000.0
wfirepmtrs.DD <- spTransform(wfirepmtrs.rpj,CRS('+init=epsg:4326'))

# Filter events based on 300 Acre cutoff
wfirepmtrs.DD <- subset(wfirepmtrs.DD, Acres >= 50.0)

# Filter errors on Oklahoma - obvious prescribed burns coded as unknown
wfirepmtrs.DD <- subset(wfirepmtrs.DD, !(STATE %in% c('OK','KS') & Fire_Type %in% c('Unknown')))
levels(droplevels(wfirepmtrs_shp$Fire_Type))
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

# Save current output, overwrite existing with idin appended
writeOGR(obj=wfirepmtrs.DD, dsn=final_path, layer='fltrdwfirepmtrs1984-2019', driver="ESRI Shapefile", overwrite=TRUE)
# Write test file to csv for use with GIS/Excel
testfile.DT <- as.data.table(wfirepmtrs.DD@data)
testfile.csv <- paste0(final_path,'/','FirePrmtrs.csv')
fwrite(testfile.DT, testfile.csv)

# Test event parameters
# Fire_ID	Fire_Name	Year	Fire_Type	Acres	STATE	idin	MARGUEEEVNTS
# CA3442911910020171205	THOMAS	              2017	Wildfire	281988	CA	13927	1
# CA3293911676620031025	CEDAR	                2003	Wildfire	268367	CA	 1061	1
# 2018-CABTU-016737	CAMP	                    2018	L6BH	    153338	CA	14676	1
# TX3015109722520110904	BASTROPCOUNTYCOMPLEX	2011	Wildfire	 31839	TX	11462	1
# CO3888410493320120623	WALDO CANYON	        2012	Wildfire	 20113	CO	 2674	1
# TN3563108347820161123	CHIMNEY TOPS 2	      2016	Wildfire	 14999	TN	11157	1
# AZ3328610958720130607	FOURMILE	            2013	Wildfire	 12961	AZ	  533	1
# CO4005110538520100906	FOURMILE CANYON	      2010	Wildfire	  5865	CO	 2723	1
# CA3859812261820171009 TUBBS                 2017  Wildfire	 36976	CA	13964	1

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

# Loop through wfirepmtrs.DD perimeter polygons
result.dt[, POLYidin := 0]
for (i in 1:nrow(wfirepmtrs.DD)) {
  # Create output directory
  POLYidin.shp <- wfirepmtrs.DD@data[i,'idin']
  result.dt[i, POLYidin := POLYidin.shp]
}
result.dt[,TEST := idin - POLYidin]
test <- result.dt[TEST != 0]
polylst <- as.vector(result.dt[COUNT > 0][, idin])

# Hand test polygons if necessary
test2 <- wfirepmtrs.DD[wfirepmtrs.DD$idin == 14676,]
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

# Hand test polygons if necessary
test2 <- wfirepmtrs.DD.fltdrd[wfirepmtrs.DD.fltdrd$idin == 14676,]
bbin <- st_bbox(test2)
xminbb = bbin[[1]]
xmaxbb = bbin[[3]]
yminbb = bbin[[2]]
ymaxbb = bbin[[4]]
xlimin = c(xminbb, xmaxbb) 
ylimin = c(yminbb, ymaxbb) 
plot(WUI.msk, xlim = xlimin, ylim = ylimin)
plot(test2, add = TRUE)

# Read attributes
wfirepmtrs.DD.fltdrd.df <- as.data.frame(wfirepmtrs.DD.fltdrd@data)

# Ensure poorly formed polygons are fixed
wfirepmtrs.DD.fltdrd.smpl <- gSimplify(wfirepmtrs.DD.fltdrd, tol = 0.00001, topologyPreserve=TRUE)

# Create a spatial polygon data frame 
wfirepmtrs.smpl <- SpatialPolygonsDataFrame(wfirepmtrs.DD.fltdrd.smpl, wfirepmtrs.DD.fltdrd.df)

# Ensure complete polygons
wfirepmtrs.DD.fltdrd.fin <- spTransform(wfirepmtrs.smpl,crs(laea_proj4))
wfirepmtrs.DD.fltdrd.buf <- gBuffer(wfirepmtrs.DD.fltdrd.fin, byid=TRUE, width=0)
wfirepmtrs.DD.fltdrd.buf.out <- spTransform(wfirepmtrs.DD.fltdrd.buf,CRS('+init=epsg:4326'))

# Test for bad polygons
sum(gIsValid(wfirepmtrs.DD.fltdrd.buf.out, byid=TRUE)==FALSE)

# Miscellaneous testing
#testout <- st_as_sfc(wfirepmtrs.DD.fltdrd.buf.out)
#testout <- lwgeom::st_make_valid(testout)
#testout.vd <- st_is_valid(testout)
#which(testout.vd == FALSE, arr.ind = TRUE)
#sp.clean <- clgeo_Clean(wfirepmtrs.DD.fltdrd) 
#report.clean <- clgeo_CollectionReport(sp.clean) 
#clgeo_SummaryReport(report.clean) 

# Hand test polygons if necessary
test2 <- wfirepmtrs.DD.fltdrd.buf.out[wfirepmtrs.DD.fltdrd.buf.out$idin == 13927,]
bbin <- st_bbox(test2)
xminbb = bbin[[1]]
xmaxbb = bbin[[3]]
yminbb = bbin[[2]]
ymaxbb = bbin[[4]]
xlimin = c(xminbb, xmaxbb) 
ylimin = c(yminbb, ymaxbb) 
plot(WUI.msk, xlim = xlimin, ylim = ylimin)
plot(test2, add = TRUE)

# Save current output
writeOGR(obj=wfirepmtrs.DD.fltdrd.buf.out, dsn=final_path, layer='final_wfirepmtrs1984-2019', driver="ESRI Shapefile", overwrite=TRUE)

# Output geojson final files post validate and fix
geojson_write(wfirepmtrs.DD.fltdrd.buf.out, lat = NULL, lon = NULL, geometry = 'polygon',
              group = NULL, file = paste0(final_path,'/','finalwfirepmtrs1984-2019_clnd.geojson'), overwrite = TRUE,
              precision = 6, convert_wgs84 = FALSE, crs = CRS('+init=epsg:4326'))

##################################################################################################################
## Section: Read accumulated exposures data aggregated through BQuery for each polygon and clean
##################################################################################################################

# Read exposures data aggregated by perimeter and 100 m grid cell
setNumericRounding(2)
getNumericRounding()
final_path <- file.path(paste0('./',DATA_DIR))
expofile <- paste0(final_path,'/','fire_perimeter_cells_total_20200408_unscaled.rds')
expo.DT <- readRDS(expofile)
setDT(expo.DT) # 13702757

#expo.DT <- expo.DT[fire_id == 'CA3859812261820171009']

# Export Thomas fire for testing
#testfile.csv <- paste0(final_path,'/','TUBBS2017_FirePrmtrs.csv')
#fwrite(expo.DT, testfile.csv)

# Re-read wild fire perimeters file exposures aggregations were based on
final_path <- file.path(paste0('./',OUTPUT_DIR))
wfirepmtrs_shp <- readOGR(paste0(final_path,'/','final_wfirepmtrs1984-2019.shp')) 
wfirepmtrs_shp <- spTransform(wfirepmtrs_shp, CRS('+init=epsg:4326'))

# Get attribute data from boundary file
wfirepmtrs_shp.DT <- as.data.table(wfirepmtrs_shp@data)

# Merge exposures data with wildfire perimeter attributes
setkey(wfirepmtrs_shp.DT, Fire_ID)
setkey(expo.DT, fire_id)
expo.DT <- merge(expo.DT, wfirepmtrs_shp.DT, by.x='fire_id', by.y='Fire_ID', all = FALSE)

# Drop obsolete columns
expo.DT[ ,`:=`(Fire_Type = NULL, Year = NULL)] 

# Sum exposures values
cols <- colnames(expo.DT)
colstosum <- cols[c(10,12)]
expo.DT[, TOTL_SM := Reduce(`+`, .SD), .SDcol = colstosum]
expo.DT <- expo.DT[TOTL_SM != 0] #2791234

# Test data Extract Camp fire
tomas.dt <- expo.DT[fire_id == 'CA3442911910020171205']

# Hand test polygons if necessary
testfire <- wfirepmtrs_shp[wfirepmtrs_shp$Fire_ID == 'CA3442911910020171205',]
bbin <- st_bbox(testfire)
xminbb = bbin[[1]]
xmaxbb = bbin[[3]]
yminbb = bbin[[2]]
ymaxbb = bbin[[4]]
xlimin = c(xminbb, xmaxbb) 
ylimin = c(yminbb, ymaxbb) 
plot(testfire, xlim = xlimin, ylim = ylimin)
points(tomas.dt[,X], tomas.dt[,Y], pch = 20, cex = 0.25, col = "green") 

# Export Thomas fire for testing
testfile.csv <- paste0(final_path,'/','THOMAS2018_FirePrmtrs.csv')
fwrite(tomas.dt, testfile.csv)

# Rename columns to be more succinct
setnames(expo.DT, 'parent_cell', 'PARNT')
setnames(expo.DT, 'cell_id', 'CELID')
setnames(expo.DT, 'geoid_county', 'CTYFPS')
setnames(expo.DT, 'exposure_subtype', 'CLOB')
setnames(expo.DT, 'building_count', 'RISKS_SM')
setnames(expo.DT, 'exposure_content_value', 'CTVAL_SM')
setnames(expo.DT, 'square_footage', 'SQFT')
setnames(expo.DT, 'exposure_value', 'BDVAL_SM')

# Create unique cell id
expo.DT <- expo.DT[,PARNT:= PARNT*1000000]
expo.DT[,EXPOID:= PARNT + CELID]
expo.DT[ ,`:=`(PARNT = NULL, CELID = NULL)]

# Summarize by event parameters
colstosum <- str_subset(names(expo.DT),'_SM')[1:4]
colstornk.indx <- match(colstosum, unlist(names(expo.DT)))
#colstogrp <- c('fire_id','year','CLOB','Fire_Name','STATE')
colstogrp <- c('fire_id')
expo.DT.evnt <- expo.DT[, lapply(.SD, sum, na.rm=TRUE), by=colstogrp, .SDcols=colstosum, drop=FALSE] 
#setorder(expo.DT.evnt,-RISKS_SM,STATE,year)
options(scipen=999)
print(expo.DT.evnt, digits=4)
#expo.DT.evnt[STATE == 'CA']

# Drop events less than $500,000
expo.DT.evnt <- expo.DT.evnt[TOTL_SM >= 500.0] 
unievnts <- unique(expo.DT.evnt, by = 'fire_id' ) # 2527 events
unievntstokeep <- as.list(unievnts[, .(fire_id)])
setkey(expo.DT, fire_id)
expo.DT <- expo.DT[unievntstokeep]

# Write exposures data.table 
final_path <- file.path(paste0('./',OUTPUT_DIR))
expofile <- paste0(final_path,'/','fire_perimeter_cells_total_cleaned.rds')
saveRDS(expo.DT, expofile)

##################################################################################################################
## Section: Investigate cleaned wildfire perimeter aggreagtions
##################################################################################################################

# Re-read combined input data
final_path <- file.path(paste0('./',OUTPUT_DIR))
expofile <- paste0(final_path,'/','fire_perimeter_cells_total_cleaned.rds')
expo.DT <- readRDS(expofile) 

# Summarize by event parameters
colstosum <- str_subset(names(expo.DT),'_SM')[1:4]
colstornk.indx <- match(colstosum, unlist(names(expo.DT)))
colstogrp <- c('fire_id','year','CLOB','Fire_Name','STATE')
expo.DT.evnt <- expo.DT[, lapply(.SD, sum, na.rm=TRUE), by=colstogrp, .SDcols=colstosum, drop=FALSE] 
setorder(expo.DT.evnt,-RISKS_SM,STATE,year)
options(scipen=999)
print(expo.DT.evnt, digits=4)
expo.DT.evnt[STATE == 'CA']

# Rank annual aggregate sums by state, clob and year - CREATE AGGREGATE LOSS DATA PER YEAR
colstornk <- str_subset(names(expo.DT.evnt),'_SM')[1:4]
colstornk.indx <- match(colstornk, unlist(names(expo.DT.evnt)))
colstogrp <- c('STATE','year','CLOB')

# Sum totals by year
annlsums <- expo.DT.evnt[, lapply(.SD, sum, na.rm=TRUE), by=colstogrp, .SDcols=colstornk] 

# calculate the rank of annual totals by state
colstogrp <- c('STATE','CLOB')
annlsums[, paste0('rank_', colstornk) := lapply(.SD, frankv, ties.method = 'min', order = -1L), 
             .SDcols = colstornk, by = colstogrp][]
setorder(annlsums,STATE,year,CLOB)
annlsums[STATE == 'CA'][CLOB == 'RES']

# Rank events by state, clob and year - CREATE OCCURENCE LOSS DATA MAX EVENT PER YEAR
#colstornk <- str_subset(names(expo.DT.evnt),'_SM')[1:4]
#colstornk.indx <- match(colstornk, unlist(names(expo.DT.evnt)))
#colstogrp <- c('STATE','year','CLOB')
#expo.DT.evnt.rnkd <- expo.DT.evnt[, .SD[frankv(.SD, ties.method="first")[.N],], by=colstogrp, .SDcols=colstornk,  drop=FALSE]
#setorder(expo.DT.evnt.rnkd,STATE,year)
#options(scipen=999)
#print(expo.DT.evnt.rnkd, digits=4)
#expo.DT.evnt.rnkd[STATE == 'CA'][CLOB == 'RES'][year == 2003]

#expo.DT.evntrnkd2 <- expo.DT.evnt[, paste0("rank_", colstornk) := lapply(.SD, frankv, ties.method = "min", order = -1L), 
#   .SDcols = colstornk, by = colstogrp][]
#expo.DT.evntrnkd2 <- expo.DT.evntrnkd2[STATE == 'CA'][CLOB == 'RES'][year == 2003]
#setorder(expo.DT.evntrnkd2,STATE,year,CLOB)
#expo.DT.evntrnkd2[, head(.SD, 1), .(STATE,year,CLOB), .SDcols=colstornk][]

# Rank events by state, clob and year - CREATE OCCURENCE LOSS DATA MAX EVENT PER YEAR
annlmaxs <- copy(expo.DT.evnt)
colstornk <- str_subset(names(annlmaxs),'_SM')[1:4]
colstornk.indx <- match(colstornk, unlist(names(annlmaxs)))
colstogrp <- c('STATE','year','CLOB')

# calculate the rank within each `id`
annlmaxs[, paste0('rank_', colstornk) := lapply(.SD, frankv, ties.method = 'min', order = -1L), 
   .SDcols = colstornk, by = colstogrp][]
annlmaxs[STATE == 'CA'][year == 2017][Fire_Name == 'TUBBS'][CLOB == 'RES']

##################################################################################################################
## Section: Filter cells by burn scar severity
##################################################################################################################

# Set data root directory
initial_path <- file.path(paste0('./',DATA_DIR))

# Initialzie CRS from raster
yrlybscr.msk <- raster::raster(paste0(initial_path,'/composite_data/MTBS_BSmosaics/','mtbs_CA_2017.tif'))    
crsin <- crs(yrlybscr.msk)
ptscrsin <- CRS('+init=epsg:4326')

# Re-read combined input data
final_path <- file.path(paste0('./',OUTPUT_DIR))
expofile <- paste0(final_path,'/','fire_perimeter_cells_total_cleaned.rds')
expo.DT <- readRDS(expofile) 

# Get the unique list of year from th einput expsours file expofile
uniyrs <- unique(expo.DT, by = 'year')
yearsin <- sort(uniyrs[,year], decreasing=FALSE)
yearsin <- yearsin[1:34]

# Test events selection
expo.DT.in <- copy(expo.DT)
#expo.DT.in <- expo.DT.in[fire_id == 'CA3859812261820171009']
#expo.DT.in <- expo.DT.in[RISKS_SM >=0.025]

# Summarize by event parameters
expo.DT.in.test <- expo.DT.in[fire_id == 'CA3859812261820171009']
colstosum <- str_subset(names(expo.DT.in.test),'_SM')[1:4]
colstornk.indx <- match(colstosum, unlist(names(expo.DT.in.test)))
colstogrp <- c('fire_id')
expo.DT.in.evnt <- expo.DT.in.test[, lapply(.SD, sum, na.rm=TRUE), by=colstogrp, .SDcols=colstosum, drop=FALSE] 
options(scipen=999)
print(expo.DT.in.evnt, digits=4)
expo.DT.in.evnt[]

# Initialize output data.table for filtered locations
expo.DT.fltrd <- data.table(NULL)

# Loop through years of exposures 
for (iyr in 1:length(yearsin)) { # iyr <- 34
  # Evaluate year from list
  iyear <- as.integer(yearsin[iyr])
  
  for (ist in 1:length(state.abb)) { # ist <- 5
  # Evaluate state from list
  state <- (state.abb[ist])
  
  # Get data for this year
  iyrdata.dt <- expo.DT.in[year == iyear][STATE == state]
  iyrdata.dt[, ID := .I] # Add id column

    # Select subset of pyromes
    if (nrow(iyrdata.dt) != 0) { 
      # prepare the 3 components: coordinates, data, and proj4string
      coords <- iyrdata.dt[ , c('X','Y')]   # coordinates
      data   <- iyrdata.dt[ , 1:dim(iyrdata.dt)[2]]          # data
      crs    <- ptscrsin # proj4string of coords
      
      # make the spatial points data frame object
      spdf <- SpatialPointsDataFrame(coords = coords,
                                     data = data, 
                                     proj4string = crs)
      spdf_trnsfrmd = spTransform(spdf,crsin)
      
      bufferedPoints <- gBuffer(spdf_trnsfrmd, width=50, byid=TRUE)
      bufferedPoints@data$idin <- attr(bufferedPoints@data, 'row.names') # Returns integers
      
      buffile <- paste0(final_path,'/bufferdata/','BufPts_',state,'_',iyear,'.shp') 
      writeOGR(obj=bufferedPoints, dsn=buffile, layer='bufferedPoints', driver='ESRI Shapefile')      
      
      # Read burn scar raster for this year
      rasfile <- paste0(initial_path,'/composite_data/MTBS_BSmosaics/','mtbs_',state,'_',iyear,'.tif')  
      
      if (file.exists(rasfile)) {
        # Read burn scar raster for this year
        yrlybscr.msk <- raster::raster(paste0(initial_path,'/composite_data/MTBS_BSmosaics/','mtbs_',state,'_',iyear,'.tif'))  
        yrlybscr.msk.CPY <- yrlybscr.msk
   
        # Create velox raster
        ras.vx <- velox(yrlybscr.msk)
    
        # Get maximum burn scar value within 150m of exposures location
        #Mode <- function(x, na.rm = FALSE) {
          #if(na.rm){
            #x = x[!is.na(x)]
          #}
          #ux <- unique(x)
          #return(ux[which.max(tabulate(match(x, ux)))])
        #}

        #result <- mclapply(seq_along(1), function(x){
        #  q <- ras.vx$crop(bufferedPoints);ras.vx$extract(bufferedPoints, fun=function(t) Mode(t,na.rm=T), small = TRUE)
        #})
    
        result <- mclapply(seq_along(1), function(x){
          q <- ras.vx$crop(bufferedPoints);ras.vx$extract(bufferedPoints, fun=NULL, small = TRUE)
        })      
        
        # Process results 
        for (i in 1:nrow(bufferedPoints)) {  
          numbers <- as.vector(unlist(result[[1]][i]))
          numbers<-numbers[!is.na(numbers)]
          getCount <- function(x) {
            u <- unique(x);
            data.frame(
             value=u,
              count=sapply(u, function(v) { length(which(x==v)) } )
            )
          }
          countout <- getCount(numbers)
          tmp <- as.vector(seq(1:6)) 
          tmp[] <- 0
          tmp[countout[,1]] <- countout[,2]
          result[[1]][i] <- sum(c(0.05,0.15,0.75,1.0,1.0,0.0)*tmp)/sum(tmp)
        }

        result.dt <- as.data.table(Reduce(rbind, unlist(result[])))
        result.dt[ , `:=`( idin = 1:.N )]
        setnames(result.dt, 'V1', 'COUNT')
    
        # Loop through wfirepmtrs.DD perimeter polygons
        result.dt[, pntidin := 0]
        for (i in 1:nrow(bufferedPoints)) {
          # Create output directory
          pntidin.shp <- bufferedPoints@data[i,'idin']
          result.dt[i, pntidin := pntidin.shp]  iyear <- 2017
        }
        result.dt[,TEST := idin - pntidin]
        test <- result.dt[TEST != 0]
        # Remove entries where many variables are NA  
        drop_rows_all_na <- function(x, pct=1) x[!rowSums(is.na(x)) >= ncol(x)*pct,]
        result.dt <- drop_rows_all_na(result.dt,0.01)
        #polylst <- as.vector(result.dt[COUNT >= 3 & COUNT !=6][, idin])
        if (iyear == 2017) {
          polylst <- as.vector(result.dt[COUNT >= 0.01][, idin])          
        } else {
          polylst <- as.vector(result.dt[COUNT >= 0.19][, idin])          
        }
      }
      
      # Write reduced buffered file
      bufferedPoints.out <- subset(bufferedPoints, (idin %in% polylst))
      buffile <- paste0(final_path,'/bufferdata/','BufPts_',state,'_',iyear,'_red.shp') 
      writeOGR(obj=bufferedPoints.out, dsn=buffile, layer='bufferedPoints.out', driver='ESRI Shapefile')
      
      # Update exposures totals
      iyrdata.dt <- iyrdata.dt[polylst]
      
      # Remove entries where many variables are NA  
      drop_rows_all_na <- function(x, pct=1) x[!rowSums(is.na(x)) >= ncol(x)*pct,]
      iyrdata.dt <- drop_rows_all_na(iyrdata.dt,0.01)
      
      # Summarize by event parameters
      colstosum <- str_subset(names(iyrdata.dt),'_SM')[1:4]
      colstornk.indx <- match(colstosum, unlist(names(iyrdata.dt)))
      colstogrp <- c('fire_id')
      iyrdata.dt.evnt <- iyrdata.dt[, lapply(.SD, sum, na.rm=TRUE), by=colstogrp, .SDcols=colstosum, drop=FALSE] 
      options(scipen=999)
      print(iyrdata.dt.evnt, digits=4)
      iyrdata.dt.evnt[]
      
      # Accrue output data.tableof filtered locations
      filesappend <- list(expo.DT.fltrd, iyrdata.dt) 
      expo.DT.fltrd <- rbindlist(filesappend, use.names=TRUE, fill=TRUE)      
    }
  }
}

# Write exposures data.table 
final_path <- file.path(paste0('./',OUTPUT_DIR))
expofile <- paste0(final_path,'/','fire_perimeter_cells_total_cleaned_filtrd.rds')
saveRDS(expo.DT.fltrd, expofile)

# Hand test polygons if necessary
testfire <- wfirepmtrs_shp[wfirepmtrs_shp$Fire_ID == 'CA3859812261820171009',]
bbin <- st_bbox(testfire)
xminbb = bbin[[1]]
xmaxbb = bbin[[3]]
yminbb = bbin[[2]]
ymaxbb = bbin[[4]]
xlimin = c(xminbb, xmaxbb) 
ylimin = c(yminbb, ymaxbb) 
plot(testfire)
plot(bufferedPoints, add = TRUE)


# Re-read combined input data
final_path <- file.path(paste0('./',DATA_DIR))
expofile <- paste0(final_path,'/','fire_perimeter_cells_damage_20200408.rds')
loss.DT <- readRDS(expofile) 
setDT(loss.DT)

##################################################################################################################
## Section: Bigquery access examples
##################################################################################################################

library(tidyverse)
library(bigrquery)
library(DBI)
library(sf)
library(mapview)
library(data.table)
options(bigrquery.page.size = 1e7)
bigrquery::bq_auth(email = "john.rowe@risq.io")
conn_rnd <- DBI::dbConnect(bigrquery::bigquery(),
                           project = "risq-rnd",
                           dataset = "test")
fire_perimeter_cells_total_20200408 <- tbl(conn_rnd, "fire_perimeter_cells_total_20200408") %>% 
  dplyr::collect()
saveRDS(fire_perimeter_cells_total_20200408, "fire_perimeter_cells_total_20200408.rds")
fire_perimeter_cells_damage_20200408 <- tbl(conn_rnd, "fire_perimeter_cells_20200408") %>% 
  dplyr::collect()
saveRDS(fire_perimeter_cells_damage_20200408, "fire_perimeter_cells_damage_20200408.rds")
fire_perimeter_cells_total_20200408 <- tbl(conn_rnd, "fire_perimeter_cells_total_20200408_unscaled") %>% 
  dplyr::collect()
saveRDS(fire_perimeter_cells_total_20200408, "fire_perimeter_cells_total_20200408_unscaled.rds")
fire_perimeter_cells_damage_20200408 <- tbl(conn_rnd, "fire_perimeter_cells_20200408_unscaled") %>% 
  dplyr::collect()
saveRDS(fire_perimeter_cells_damage_20200408, "fire_perimeter_cells_damage_20200408_unscaled.rds")

library(tidycensus)
library(tidyverse)

# Request api key http://api.census.gov/data/key_signup.html.
#census_api_key('92e89f1664544e58c6c6c3f4d289a5e9455d4e4b', install=TRUE)
census_api_key('92e89f1664544e58c6c6c3f4d289a5e9455d4e4b'), overwrite=TRUE)
readRenviron("~/.Renviron")

v17 <- load_variables(2017, "acs5", cache = TRUE)
write.table(v17 , file = 'C:/GitRepositories/WFireOutcomeWork/ACSAttrbts.csv')

MedHouseDolrs <- get_acs(geography = 'county', 
                         variables = c(medianhseprce = 'B25077_001'),  
                         year = 2018)
write.table(MedHouseDolrs , file = 'C:/GitRepositories/WFireOutcomeWork/ACSAttrbts_HsePrice.csv')

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








