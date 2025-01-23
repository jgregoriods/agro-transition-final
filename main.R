required_packages <- c(
    "ape", "caret", "dplyr", "fastshap", "fields", "ggplot2", "gridExtra", 
    "latticeExtra", "progress", "randomForestSRC", "rasterVis", "RColorBrewer",
    "rnaturalearth", "rnaturalearthdata", "sf", "sp", "shapviz", "terra", 
    "tidyterra", "vegan", "ggpubr", "ade4"
)

missing_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(missing_packages)) {
    print("Installing dependencies...")
    install.packages(missing_packages)
}

library(ape)
library(caret)
library(dplyr)
library(fastshap)
library(fields)
library(ggplot2)
library(gridExtra)
library(latticeExtra)
library(progress)
library(randomForestSRC)
library(rasterVis)
library(RColorBrewer)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(sp)
library(shapviz)
library(terra)
library(tidyterra)
library(vegan)
library(ggpubr)
library(ade4)  # Dependency for quickMEM

sf_use_s2(FALSE)

# please download the quickMEM source code from
# https://github.com/ajsmit/Quantitative_Ecology/blob/main/Num_Ecol_R_book_ed1/quickMEM.R
source("quickMEM.R")

source("var_names.R")

# ---------------------------------------------------------
# Data preparation

# Define datum
WGS84 <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'

# Mask the NA cells for plotting
bio <- rast("rasters/bioclim/bio01.asc")
alt <- rast("rasters/terrain/altitude.asc")
land_mask <- (!is.na(bio) & !is.na(alt))

# For the figures
coast <- st_union(ne_countries())

# Read data
raw_data <- read.csv("data/c14dates.csv")

all_data_points <- st_as_sf(raw_data, coords=c("Longitude", "Latitude"), crs=WGS84)
all_data_points <- all_data_points[order(all_data_points$AgeCalBP), ]

# Ages
age_plot <- ggplot() +
    geom_sf(data=coast, fill="white") +
    geom_sf(data=all_data_points, aes(col=AgeCalBP), size=0.5, alpha=0.5) +
    scale_color_viridis_c(option="turbo", name="Age (cal BP)") +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    coord_sf(ylim = c(-60, 90))
jpeg("figs/C14Ages.jpg", width=2000, height=1000, res=300)
plot(age_plot)
dev.off()

# Sources
all_data_points$Database[all_data_points$Database == ""] <- "Other source"
n <- length(unique(all_data_points$Database))
source_plot <- ggplot() +
    geom_sf(data=coast, fill="white") +
    geom_sf(data=all_data_points, aes(col=Database), size=0.5, alpha=0.5) +
    scale_color_manual(values = colorRampPalette(brewer.pal(8, "Set3"))(n)) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    coord_sf(ylim = c(-60, 90))
jpeg("figs/C14Sources.jpg", width=2000, height=1000, res=300)
plot(source_plot)
dev.off()

# Read all environmental rasters
# we read them separately for varpart
bioclim_files <- list.files("rasters/bioclim", full.names=T)
soil_files <- list.files("rasters/edaphic", full.names=T)
terrain_files <- list.files("rasters/terrain", full.names=T)

layers <- c(rast(bioclim_files), rast(soil_files), rast(terrain_files))

# Fill NAs within a radius of 3 cells to avoid loosing data
# especially on coastal regions.
layers <- focal(layers, w=matrix(1,3,3), fun=mean, na.policy="only")

# Create layers with latitude and longitude to be added as predictors
rlon <- rlat <- rast(layers[[1]])
xy <- xyFromCell(layers[[1]], 1:ncell(layers[[1]]))
values(rlon) <- xy[,1]
values(rlat) <- xy[,2]

names(rlon) <- "xcoord"
names(rlat) <- "ycoord"

layers <- c(layers, rlon, rlat)

# Extract environmental data for each point and create a data frame
data <- terra::extract(layers, raw_data[c("Longitude", "Latitude")], ID=FALSE)

# Target variable is time since first transition
data$y <- max(raw_data$AgeCalBP) - raw_data$AgeCalBP
data <- na.omit(data)

# Store latitude and longitude for use later
xcoord <- data$xcoord
ycoord <- data$ycoord

# Predictors and target
X <- data[,1:(ncol(data)-3)]
y <- data$y

# Remove collinearity
cor_vars <- findCorrelation(cor(X), names=T)
X <- X[,!colnames(X) %in% cor_vars]
findLinearCombos(X) # none
layers <- subset(layers, cor_vars, negate=T)

data <- data[,!colnames(data) %in% cor_vars]

write.csv(data, "saved_data/data.csv", row.names=FALSE)
# data <- read.csv("saved_data/data.csv")

# ---------------------------------------------------------
# distance-based Moran's eigenvector maps (dbMEM)

jpeg("figs/SI_mem.jpg", width=2000, height=1600, res=300)
sink("results/mem.txt")
mem <- quickMEM(data$y, data[,c("xcoord", "ycoord")])
sink()
dev.off()

# Variation partitioning
terrain <- X[,paste(colnames(X), ".asc", sep="") %in% basename(terrain_files)]
bioclim <- X[,paste(colnames(X), ".asc", sep="") %in% basename(bioclim_files)]
soil <- X[,paste(colnames(X), ".asc", sep="") %in% basename(soil_files)]

spatial <- data[,c("xcoord", "ycoord")]

v <- varpart(y, terrain, bioclim, soil, spatial)

jpeg("figs/SI_varpart.jpg", width=2000, height=1600, res=300)
plot(v, Xnames=c("Terrain", "Bioclim", "Soil", "Spatial"), bg=1:4)
dev.off()

# For random forest
data$event <- 1

# ---------------------------------------------------------
# Random Forest

median_survival <- function(row) {
    idx <- which(row < 0.5)
    if (length(idx) > 0) {
        return(rf_model$time.interest[idx[1]])
    } else {
        return(NA)
    }
}

set.seed(100)
print("Tuning random forest...")
o <- tune(Surv(y, event) ~ ., data=data, sampsize=nrow(data) * 0.8)

print("Fitting random forest...")
rf_model <- rfsrc(Surv(y, event) ~ ., data=data, importance=TRUE,
                  nodesize=o$optimal[1], mtry=o$optimal[2], save.memory = TRUE,
                  sampsize=nrow(data) * 0.8)

sink("results/rsf.txt")
print(rf_model)
sink()

names(rf_model$importance) <- var_names[names(rf_model$importance)]
jpeg("figs/SI_RF.jpg", width=3000, height=2000, res=300)
plot(rf_model)
dev.off()

predicted <- apply(rf_model$survival, 1, median_survival)
residuals <- y - predicted

data_points <- st_as_sf(data, coords=c("xcoord", "ycoord"), crs=WGS84)
data_points["pred"] <- max(raw_data$AgeCalBP) - predicted
data_points["res"] <- residuals
data_points <- data_points[order(abs(data_points$res)), ]

# Individual predictions
pred_plot <- ggplot() +
    geom_sf(data=coast, fill="white") +
    geom_sf(data=data_points, aes(col=pred), size=0.5, alpha=0.5) +
    scale_color_viridis_c(option="turbo", name="Predicted yr BP") +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    coord_sf(ylim = c(-60, 90))
jpeg("figs/SI_pred.jpg", width=2000, height=1000, res=300)
plot(pred_plot)
dev.off()

# Residuals
res_plot <- ggplot() +
    geom_sf(data=coast, fill="white") +
    geom_sf(data=data_points, aes(col=res)) +
    scale_color_viridis_c(name="Residual", option="plasma") +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    coord_sf(ylim = c(-60, 90))
png("figs/SI_res.png", width=2000, height=1000, res=300)
plot(res_plot)
dev.off()

weights = rdist(data[,c("xcoord", "ycoord")])

sink("results/moran.txt")
print(Moran.I(x = residuals, w = weights))
sink()

# Convert raster stack to data frame for prediction
df <- as.data.frame(layers)
df[is.na(df)] <- 0

print("Predicting random forest...")
raw_rf_pred <- predict(rf_model, df, na.action="na.omit")

times <- c(4000, 6000, 8000, 10000)
indices <- lapply(times, function(t) {
    which.min(abs(rf_model$time.interest - t))
})

rf_stack <- NULL
for (i in indices) {
    # 1 - probability of "survival" = probability that the event happened
    time_slice <- rast(layers[[1]])
    values(time_slice) <- 1 - raw_rf_pred$survival[,i]
    if (is.null(rf_stack)) {
        rf_stack <- time_slice
    } else {
        rf_stack <- c(rf_stack, time_slice)
    }
}
names(rf_stack) <- c("8k", "6k", "4k", "2k")

# Mask by land
rf_stack[land_mask == 0] <- NA

predicted_time <- apply(raw_rf_pred$survival, 1, median_survival)
predicted_time_r <- layers[[1]]
values(predicted_time_r) <- max(raw_data$AgeCalBP) - predicted_time
predicted_time_r[land_mask == 0] <- NA

rm(raw_rf_pred)
gc()

save(rf_model, file="saved_data/rf_model.RData")
# load("saved_data/rf_model.RData")

writeRaster(rf_stack, "saved_data/rf_stack.tif", overwrite=TRUE)
# rf_stack <- rast("saved_data/rf_stack.tif")

writeRaster(predicted_time_r, "saved_data/predicted_time_r.tif", overwrite=TRUE)
# predicted_time_r <- rast("saved_data/predicted_time_r.tif")

# Random forest plots
rf_plot <- ggplot() +
    geom_sf(data=coast, fill="white") +
    geom_spatraster(data=rf_stack) +
    facet_wrap(~lyr) +
    geom_sf(data=coast, fill="transparent") +
    scale_fill_viridis_c(option="turbo", na.value="transparent", name="Probability") +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    coord_sf(ylim = c(-60, 90))
jpeg("figs/RForest.jpg", width=2000, height=1000, res=300)
plot(rf_plot)
dev.off()

# Predicted time
rf_pred <- ggplot() +
    geom_sf(data=coast, fill="white") +
    geom_spatraster(data=predicted_time_r) +
    geom_sf(data=coast, fill="transparent") +
    scale_fill_viridis_c(option="turbo", na.value="transparent", name="Predicted yr BP") +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    coord_sf(ylim = c(-60, 90))
jpeg("figs/RForestPred.jpg", width=2000, height=1000, res=300)
plot(rf_pred)
dev.off()

# ---------------------------------------------------------
# Shapley values

pfun <- function(object, newdata) {
    predict(object, newdata)$predicted
}

xvars <- select(data, -c(y, event))

print("Calculating Shapley values")
print("This may take a while...")
system.time({explainer <- explain(rf_model, X=xvars, pred_wrapper=pfun, nsim=50, adjust=TRUE)})

colnames(xvars) <- var_names[colnames(xvars)]
colnames(explainer) <- colnames(xvars)

save(explainer, file="saved_data/explainer.RData")
# load("saved_data/explainer.RData")

shv <- shapviz(explainer, X=xvars)

jpeg("figs/ShapSummary.jpg", width=2000, height=1000, res=300)
sv_importance(shv, max_display=10, kind="beeswarm", size=0.5, alpha=0.5) +
    scale_color_gradient(low="blue", high="red")
dev.off()

top_vars <- names(sort(apply(abs(explainer), 2, mean), decreasing=T))
top_vars <- top_vars[top_vars != "Longitude" & top_vars != "Latitude"]

shap_scores_df <- as.data.frame(explainer)

shap_scores_df$xcoord <- data$xcoord
shap_scores_df$ycoord <- data$ycoord
shap_points <- st_as_sf(shap_scores_df, coords=c("xcoord", "ycoord"), crs=WGS84)

# Spatial distribution of Shapley values
shap_maps <- lapply(1:4, function(i) {
    shap_points_ordered <- shap_points %>% arrange(get(top_vars[i]))
    ggplot() +
        geom_sf(data=coast, fill="white") +
        geom_sf(data=shap_points_ordered, aes(color=get(top_vars[i])), size=0.5, alpha=0.5) +
        scale_color_viridis_c(option="turbo", name="Shapley value") +
        scale_x_continuous(expand = c(0,0)) +
        scale_y_continuous(expand = c(0,0)) +
        coord_sf(ylim = c(-60, 90)) +
        labs(title=top_vars[i]) +
        theme(
            plot.title = element_text(
                size = 10,
                hjust = 0.5
            ),
            legend.title = element_text(size=8)
        )
})

jpeg("figs/ShapMaps.jpg", width=2000, height=1000, res=300)
# do.call(grid.arrange, shap_maps)
ggarrange(plotlist=shap_maps, common.legend=TRUE, legend="right")
dev.off()

dependence_plots <- lapply(1:4, function(i) {
    sv_dependence(shv, v=top_vars[i], size=0.5, alpha=0.5) +
        scale_color_gradient(low="blue", high="red") +
        # change y labels to "Shapley value"
        labs(y="Shapley value")
})

# Dependence plots
jpeg("figs/ShapDependence.jpg", width=2500, height=1800, res=300)
# sv_dependence(shv, v=top_vars[1:4])
ggarrange(plotlist=dependence_plots, common.legend=TRUE, legend="right")
dev.off()
