#################
### libraries ###
#################

library("raster")
library("tmap")
library("OpenStreetMap")
library("colorspace")
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
	### loading data ###
	metricsNames <- get_metricsNames(DirFFM)
	GagesUSGS <- read_GagesUSGS(FileGages)
	### spatial join ###
	GagesUSGS_joined <- make_spatial_join(GagesUSGS, DirFFM, metricsNames)
	### write results  ###
	if (!dir.exists("output")) dir.create("output")
	if (!dir.exists("output/shp")) dir.create("output/shp")
	shapefile(GagesUSGS_joined, file = file.path("output/shp/GagesUSGS_joined.shp"), overwrite = TRUE)
	return(GagesUSGS_joined)
}

make_maps <- function(GagesUSGS_joined, start = 2){
	get.aes <- function(name){
		if (grepl("percent_IQR", name)){
			brk <- c(0, 25, 35, 45, 55, 65, 75, 100)
			pal1 <- head(sequential_hcl(7, "Blues 2"), 3)
			pal2 <- sequential_hcl(2, "Greens 2")[1]
			pal3 <- tail(sequential_hcl(7, "Reds 2", rev = TRUE), 3)
			pal <- c(pal1, pal2, pal3)
			labs <- NULL
		} else if (grepl("pred_med_IQR", name)){
			pal1 <- head(sequential_hcl(7, "Blues 2"), 1)
			pal3 <- tail(sequential_hcl(7, "Reds 2", rev = TRUE), 1)
			pal <- c(pal3, pal1)
			brk <- NULL
			labs <- c("False (0)", "True (1)")
		} else if (grepl("medOE", name)){
			brk <- c(-Inf, 0.5, 0.7, 0.9, 1.1, 1.4, 2, Inf)
			pal1 <- head(sequential_hcl(7, "Blues 2"), 3)
			pal2 <- sequential_hcl(2, "Greens 2")[1]
			pal3 <- tail(sequential_hcl(7, "Reds 2", rev = TRUE), 3)
			pal <- c(pal1, pal2, pal3)
			labs <- NULL
		}
		return(list(pal = pal, brk = brk, labs = labs))
	}
	if (file.exists(file.path("output/map.Rds"))){
		overwriteBool <- readline(prompt="Output files already exists. Do you want to overwrite it? [TRUE/FALSE]")
	} else {
		overwriteBool <- TRUE
	}
	if (overwriteBool){
		if (!dir.exists("output/pdfs")) dir.create("output/pdfs")
		boundingBox <- tmaptools::geocode_OSM("California")$bbox
		upperLeft <- c(boundingBox$ymax, boundingBox$xmin)
		lowerRight <- c(boundingBox$ymin, boundingBox$xmax)
		osm <- raster(openmap(upperLeft, lowerRight, type = "http://tile.stamen.com/toner-background/{z}/{x}/{y}.png", zoom = 6))
		background <- tm_shape(osm) + tm_rgb(alpha = 1/3) + tm_scale_bar(position = c("right", "bottom")) + tm_compass(position = c("left", "bottom")) + tm_layout(legend.position = c("right", "top"))
		name <- names(GagesUSGS_joined)[start]
		ma <- tm_shape(GagesUSGS_joined, name = name) 
		ma <- ma + tm_dots(col = name, 
				title = name, 
				palette = get.aes(name)$pal,
				breaks = get.aes(name)$brk, 
				style = ifelse(is.null(get.aes(name)$brk), "cat", "fixed"),
				labels = get.aes(name)$labs,
				size = .25) 
		maPdf <- background + ma
		tmap_save(maPdf, dpi = 300, filename = file.path(paste0("output/pdfs/", name, ".pdf")), units = "mm", width = 297, height = 210)
		end <- ncol(GagesUSGS_joined)
		for (i in seq((start+1),end)){
			name <- names(GagesUSGS_joined)[i]
			maPdf <- tm_shape(GagesUSGS_joined, name = name)
			maPdf <- maPdf + tm_dots(col = name, 
				title = name, 
				palette = get.aes(name)$pal,
				breaks = get.aes(name)$brk, 
				style = ifelse(is.null(get.aes(name)$brk), "cat", "fixed"),
				labels = get.aes(name)$labs,
				size = .25)  
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