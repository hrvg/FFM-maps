#################
### libraries ###
#################

library("raster")
library("tmap")
library("leaflet")

##################
### parameters ###
##################

FileGages <- "data/Final_Reference_Gages/Final_Reference_Gages_7.3.2018.csv"
DirFFM <- "data/FFM Site Performance"

#################
### functions ###
#################

get_STAIDs <- function(DirFFM){
	listFiles <- list.files(file.path(DirFFM), full.names = TRUE)
	STAIDs <- sapply(listFiles, function(File){
		data <- read.csv(File)
		data$STAIDT
	})
	STAIDs <- unique(unname(unlist(STAIDs)))
	STAIDs <- gsub("T", "", STAIDs)
	return(STAIDs)
}

get_metricsNames <- function(DirFFM){
	metricsNames <- list.files(file.path(DirFFM), full.names = FALSE)
	metricsNames <- metricsNames[c(TRUE, FALSE)]
	metricsNames <- unname(sapply(metricsNames, function(name) unlist(strsplit(name, "_site"))[1]))
	return(metricsNames)
}

get_points <- function(raw_mat){
	names <- raw_mat$Name
	lon <- raw_mat$long
	lat <- raw_mat$lat
	lonlat <- cbind(lon,lat)
	crdref <- CRS("+proj=longlat +datum=WGS84")
	pts <- SpatialPoints(lonlat, proj4string = crdref)
	return(pts)
}

read_GagesUSGS <- function(FileGages){
	STAIDs <- get_STAIDs(DirFFM)
	if (tools::file_ext(FileGages) == "shp"){
		GagesUSGS <- shapefile(file.path(FileGages))
	} else {
		GagesUSGScsv <- read.csv(FileGages)
		keeps <- c("USGS_GAGE", "LATITUDE", "LONGITUDE")
		GagesUSGScsv <- GagesUSGScsv[, colnames(GagesUSGScsv) %in% keeps]
		GagesUSGScsv <- na.omit(GagesUSGScsv)
		lonlat <- cbind(GagesUSGScsv$LONGITUDE, GagesUSGScsv$LATITUDE)
		crsRef <- CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")
		GagesUSGS <- SpatialPoints(lonlat, proj4string = crsRef)
		GagesUSGS$STAID <- GagesUSGScsv$USGS_GAGE
	}
	if(table(STAIDs %in% GagesUSGS$STAID)["TRUE"] != length(STAIDs)){
		warning(paste0(table(STAIDs %in% GagesUSGS$STAID)["FALSE"], " STAID were not found. Proceeding."))
	}
	GagesUSGS <- GagesUSGS[which(GagesUSGS$STAID %in% STAIDs), ]
	return(GagesUSGS)
}

make_spatial_join <- function(GagesUSGS, DirFFM, metricsNames, performanceMetrics = c("percent_IQR", "pred_med_IQR", "medOE")){
	listFiles <- list.files(file.path(DirFFM), full.names = TRUE)
	for (metric in metricsNames){
		Files <- listFiles[grep(metric, listFiles)]
		data <- do.call(cbind, lapply(Files, read.csv))
		performanceData <- performanceMetrics[performanceMetrics %in% names(data)]  		
		if(length(performanceData) != length(performanceMetrics)){
			for (elmt in setdiff(performanceMetrics, performanceData)) warning(paste0(elmt, " not found for: ", metric,". Proceeding."))
		}
		resultsFFM <- as.data.frame(matrix(NA, nrow = length(GagesUSGS), ncol = length(performanceData)))
		for (i in seq_along(performanceData)){
			ind <- na.omit(match(gsub("T", "", data$STAIDT), GagesUSGS$STAID))
			resultsFFM[ind, i] <- data[[performanceData[i]]][which(gsub("T", "", data$STAIDT) %in% GagesUSGS$STAID)]
		}
		colnames(resultsFFM) <- paste(metric, performanceData, sep = ".")
		GagesUSGS@data <- data.frame(GagesUSGS@data, resultsFFM)
	}
	return(GagesUSGS)
}

get_GagesUSGS_joined <- function(DirFFM, FileGages){
	if (file.exists(file.path("output/GagesUSGS_joined.shp"))){
		overwriteBool <- readline(prompt="Joined shapefile already exists. Do you want to overwrite it? [TRUE/FALSE]")
	} else {
		overwriteBool <- TRUE
	}
	if (overwriteBool){
		### loading data ###
		metricsNames <- get_metricsNames(DirFFM)
		GagesUSGS <- read_GagesUSGS(FileGages)
		### spatial join ###
		GagesUSGS_joined <- make_spatial_join(GagesUSGS, DirFFM, metricsNames)
		### write results  ###
		if (!dir.exists("output")) dir.create("output")
		shapefile(GagesUSGS_joined, file = file.path("output/GagesUSGS_joined.shp"), overwrite = TRUE)
	} else {
		GagesUSGS_joined <- shapefile(file.path("output/GagesUSGS_joined.shp"))
	}
	return(GagesUSGS_joined)
}

make_maps <- function(GagesUSGS_joined, start = 15, end = 70){
	if (file.exists(file.path("output/map.Rds"))){
		overwriteBool <- readline(prompt="Output files already exists. Do you want to overwrite it? [TRUE/FALSE]")
	} else {
		overwriteBool <- TRUE
	}
	if (overwriteBool){
		ma <- tm_shape(GagesUSGS_joined, name = names(GagesUSGS_joined)[start]) 
		ma <- ma + tm_dots(col = names(GagesUSGS_joined)[start], title = names(GagesUSGS_joined)[start], palette = "-viridis") 
		for (i in seq((start+1),end)){
			name <- names(GagesUSGS_joined)[i]
			ma <- ma + tm_shape(GagesUSGS_joined, name = name)
			ma <- ma + tm_dots(col = name, 
				title = name, 
				palette = ifelse(length(unique(GagesUSGS_joined[[name]])) < 3, "Set1", "-viridis"))
		}
		ma <- ma + tm_basemap(server = "Esri.WorldTopoMap")
		saveRDS(names(GagesUSGS_joined)[start:end], file.path("output/layer_names.Rds"))
		saveRDS(ma, file.path("output/map.Rds"))
	} else {
		ma <- readRDS(file.path("output/map.Rds"))
	}	
	return(ma)
}

####################
### spatial join ###
####################

GagesUSGS_joined <- get_GagesUSGS_joined(DirFFM, FileGages)

############
### maps ###
############

if (tools::file_ext(FileGages) == "shp"){
	ma <- make_maps(GagesUSGS_joined)
} else {
	ma <- make_maps(GagesUSGS_joined, start = 2, end = 56)
}