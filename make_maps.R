#################
### libraries ###
#################

library("raster")
library("tmap")

##################
### parameters ###
##################

FileGages <- "data/gagesII_9322_point_shapefile/gagesII_9322_sept30_2011.shp"
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

read_GagesUSGS <- function(FileGages){
	GagesUSGS <- shapefile(file.path(FileGages))
	STAIDs <- get_STAIDs(DirFFM)
	if(table(STAIDs %in% GagesUSGS$STAID)["FALSE"] != 0){
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

####################
### spatial join ###
####################

GagesUSGS_joined <- get_GagesUSGS_joined(DirFFM, FileGages)

############
### maps ###
############

listMaps <- lapply(seq(16, 70), function(i){
	name <- names(GagesUSGS_joined)[i]
	print(name)
	ma <- tm_shape(GagesUSGS_joined, name = name) 
	ma <- ma + tm_dots(col = name, title = name, palette = "-viridis")
	ma <- ma + tm_basemap(server = "Esri.WorldTopoMap")
	return(ma)
})