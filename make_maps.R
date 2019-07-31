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
		crsRef <- CRS("+proj=longlat +datum=WGS84")
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
	if (file.exists(file.path("output/shp/GagesUSGS_joined.shp"))){
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
		if (!dir.exists("output/shp")) dir.create("output/shp")
		shapefile(GagesUSGS_joined, file = file.path("output/shp/GagesUSGS_joined.shp"), overwrite = TRUE)
	} else {
		GagesUSGS_joined <- shapefile(file.path("output/shp/GagesUSGS_joined.shp"))
	}
	return(GagesUSGS_joined)
}

make_maps <- function(GagesUSGS_joined, start = 2){
	if (file.exists(file.path("output/map.Rds"))){
		overwriteBool <- readline(prompt="Output files already exists. Do you want to overwrite it? [TRUE/FALSE]")
	} else {
		overwriteBool <- TRUE
	}
	if (overwriteBool){
		if (!dir.exists("output/pdfs")) dir.create("output/pdfs")
		osm <- tmaptools::read_osm(tmaptools::bb(GagesUSGS_joined))
		background <- tm_shape(osm) + tm_rgb() + tm_scale_bar(position = c("right", "bottom")) + tm_compass(position = c("left", "bottom")) + tm_layout(legend.position = c("right", "top"))
		name <- names(GagesUSGS_joined)[start]
		ma <- tm_shape(GagesUSGS_joined, name = name) 
		ma <- ma + tm_dots(col = name, 
				title = name, 
				palette = ifelse(length(unique(GagesUSGS_joined[[name]])) < 3, "Set1", "-viridis"),
				size = .25,
				alpha = .8) 
		maPdf <- background + ma
		tmap_save(maPdf, dpi = 300, filename = file.path(paste0("output/pdfs/", name, ".pdf")), units = "mm", width = 297, height = 210)
		end <- ncol(GagesUSGS_joined)
		for (i in seq((start+1),end)){
			name <- names(GagesUSGS_joined)[i]
			maPdf <- tm_shape(GagesUSGS_joined, name = name)
			maPdf <- maPdf + tm_dots(col = name, 
				title = name, 
				palette = ifelse(length(unique(GagesUSGS_joined[[name]])) < 3, "Set1", "-viridis"),
				size = .25,
				alpha = .8) 
			ma <- ma + maPdf
			maPdf <- background + maPdf
			tmap_save(maPdf, dpi = 300, filename = file.path(paste0("output/pdfs/", name, ".pdf")), units = "mm", width = 297, height = 210)
		}
		ma <- ma 
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
	ma <- make_maps(GagesUSGS_joined, start = 15)
} else {
	ma <- make_maps(GagesUSGS_joined, start = 2)
}

# osmSFE <- tmaptools::read_osm(rrSFE)
# ma1 <- tm_shape(osmSAC) + tm_rgb() + tm_shape(rrSAC, name = "NC uncertainty hotspots") + tm_raster(palette = "-inferno", n = ncol, contrast = c(.2, .8), alpha = .7, breaks = brk, legend.show = FALSE, title = "Intensity") + tm_scale_bar(position = c("right", "bottom")) + tm_compass(position = c("left", "bottom")) 
# ma2 <- tm_shape(osmSFE) + tm_rgb() + tm_shape(rrSFE, name = "NC uncertainty hotspots") + tm_raster(palette = "-inferno", n = ncol, contrast = c(.2, .8), alpha = .7, breaks = brk, legend.show = TRUE, legend.format = list(format = "f", digits = 1), title = "Intensity") + tm_scale_bar(position = c("right", "top")) + tm_compass(position = c("right", "top")) + tm_layout(legend.position = c("left", "bottom")) + tm_basemap(server = "Esri.WorldTopoMap")
# tm <- tmap_arrange(ma1, ma2)
# tmap_save(ma1, dpi = 900, filename = "SACmap.pdf", units = "mm", width = 297, height = 210)