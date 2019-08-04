# Requirements

This code ran with `R` v3.6

The following packages are required for generating the output:

- `raster`
- `tmap`
- `OpenStreetMap`
- `colorspace`

The following packages are required for visualizing the output with a shinyApp:

- `leaflet`
- `shiny`
- `flexdashboard`

# Scripts

- `make_maps.R`: generates the output
- `FFM_maps.Rmd`: is a lightweight shinyApp to visualize the results (not pushed online)

# Outputs

- the `.Rds` files contain the `tmap` objects
- `output/shp` contains the shapefile with the gages coordinates associated with the available metrics at each gage
- `output/pdfs` contains `.pdf` maps for each metric