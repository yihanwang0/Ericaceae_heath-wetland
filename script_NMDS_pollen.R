library(vegan)
library(cluster) # for Gower distance
library(ggplot2) # for custom plotting
library(dplyr)
library(ggrepel) # automatically adjust label positions to avoid overlap
# install.packages("ggrepel")

#load dataset
E.data <- read.csv("dataset2_for_script.csv", header = TRUE, fileEncoding = "ISO-8859-1") # header = TRUE : CSV contains cloumn names (headers)
# colnames(E.data)

# Filter fossil data
Fossil <- E.data %>% filter(Type != "Reference") # fossil data only
Reference <- E.data %>% filter(Type == "Reference") # reference data only

# treat coded qualitative variables as factors
E.data$`Dispersal.unit` <- as.factor(E.data$`Dispersal.unit`)
E.data$`Isodiametric.tetrad.or.anisodiametric.tetrad` <- as.factor(E.data$`Isodiametric.tetrad.or.anisodiametric.tetrad`)
E.data$`Aperture.type` <- as.factor(E.data$`Aperture.type`)
E.data$`Amb.shape..polar.view.shape.` <- as.factor(E.data$`Amb.shape..polar.view.shape.`)
E.data$`Exine.ornamentation..surface.` <- as.factor(E.data$`Exine.ornamentation..surface.`)
E.data$`Equatorial.shape` <- as.factor(E.data$`Equatorial.shape`)

# select the numerical and factor variables for NMDS
quantitatives <- E.data %>% select(where(is.numeric))
qualitatives <- E.data %>% select (where(is.factor))

# combine numerics and factors for Gower distance calculation
combined <- bind_cols(quantitatives, qualitatives)

# calculate Gower distance for this mix dataset
gower_distance <- daisy(combined, metric = "gower")

# perform nmds using gower_distance
nmds <- metaMDS(gower_distance, k = 2, trymax = 250)

# extract nmds coordinates (scores) for each grain, convert nmds scores into a data frame
nmds_coordinates <- as.data.frame(scores(nmds))

# Save NMDS coordinates to a CSV file
write.csv(nmds_coordinates, "nmds_coordinates.csv", row.names = FALSE)

# Load the NMDS coordinates from the CSV file
nmds_coordinates <- read.csv("nmds_coordinates.csv")

# add species, genus, fossil type info as additional columns to the data frame for labeling
nmds_coordinates$`Genus` <- E.data$`Genus`
nmds_coordinates$`Species` <- E.data$`Species`
nmds_coordinates$`Fossil.type` <- E.data$`Fossil.type`

# Save a CSV file
# write.csv(nmds_coordinates, "nmds_coordinates_additional_columns_for_labeling.csv", row.names = FALSE)
# Load the NMDS coordinates from the CSV file
# nmds_coordinates <- read.csv("mds_coordinates_additional_columns_for_labeling.csv")

# Separate fossil coordinates (do not include genus/species for fossils; they will be NA)
fossil_coordinates <- nmds_coordinates %>% 
  filter(is.na(Genus) & is.na(Species))  # Select rows where Genus and Species are NA (fossil data)

# Separate reference coordinates (to ensure fossil data does not interfere with convex hulls)
reference_coordinates <- nmds_coordinates %>% 
  filter(!is.na(Genus) & !is.na(Species))  # Select rows where Genus and Species are not NA (reference data)

# extract stress value, round to 3 decimal points (below 0.2 indicate a good fit)
stress <- round (nmds$stress, 4)

# Calculate centroids (the average positions) of the points for each species
centroids_species <- reference_coordinates %>%
  group_by(Species) %>%
  summarise(NMDS1 = mean(NMDS1), NMDS2 = mean(NMDS2))

# Calculate centroids (the average positions) of the points for each genus
centroids_genus <- reference_coordinates %>%
  group_by(Genus) %>%
  summarise(NMDS1 = mean(NMDS1), NMDS2 = mean(NMDS2))

# Calculate genus-level convex hulls (polygons for genus clusters)
genus_hulls <- reference_coordinates %>%
  group_by(Genus) %>%
  slice(chull(NMDS1, NMDS2))  # Select points forming the convex hull for each genus

# Calculate species-level convex hulls (polygons for species clusters)
species_hulls <- reference_coordinates %>%
  group_by(Species) %>%
  slice(chull(NMDS1, NMDS2))  # Select points forming the convex hull for each species

# --------------------------------------

# Plot NMDS clusters
ggplot(nmds_coordinates, aes(x = NMDS1, y = NMDS2)) +

# Plot genus-level filled convex hull
geom_polygon(data = genus_hulls, aes(x = NMDS1, y = NMDS2, fill = Genus), alpha = 0.2, color = NA) +
  
# Plot species-level convex hull (outlined shape, no fill)
geom_polygon(data = species_hulls, aes(x = NMDS1, y = NMDS2, group = Species, color = Species), fill = NA, linetype = "dashed") +
 
# Add fossil data as red dots (with no genus/species labels)
geom_point(data = fossil_coordinates, aes(x = NMDS1, y = NMDS2, shape = Fossil.type), size = 3, color = "red") +

# Customize shape for different fossil types
scale_shape_manual(values = c(0, 1, 2)) + # Assign different symbols for each type

# Add genus labels at centroids with automatic repelling to avoid overlap (only one label per genus)
geom_text_repel(data = centroids_genus, aes(x = NMDS1, y = NMDS2, label = Genus), size = 6, hjust = 0, vjust = 1.5, fontface = "bold", nudge_y = -0.007) +

# Add species labels at centroids with automatic repelling to avoid overlap (only one label per species)
geom_text_repel(data = centroids_species, aes(x = NMDS1, y = NMDS2, label = Species), size = 4, hjust = 0, vjust = 1.5, fontface = "bold", max.overlaps = 1) +
  
# Add stress value to the title
labs(title = paste("NMDS Ordination (Gower Distance) of Reference and Fossil Pollen Morphometrics ", "\nStress =", stress),
       x = "NMDS1", y = "NMDS2", shape = "Fossil Type") +

# for fill colors: set the color palette for the filled polygons that represent genus clusters
scale_fill_manual(values = rainbow(length(unique(genus_hulls$Genus))))+

# use a color palette for species outlines, assign a unique color to each species's outline
scale_color_manual(values=rainbow(length(unique(species_hulls$Species)))) +

# Theme clear for output
theme_classic() +

  # Customize theme to show black axis lines (x and y axes)
  theme(
    # Add axis lines (black lines along x and y axes)
    axis.line.x = element_line(color = "black", size = 1),  # Black x-axis
    axis.line.y = element_line(color = "black", size = 1),  # Black y-axis
    
    # keep axis ticks
    axis.ticks = element_line(color = "black", size = 1),
    axis.ticks.length = unit(0.5, "cm")  # Customize tick length
    )+
    
theme(
    axis.text.x = element_text(size = 12),  # Adjust x-axis text size
    axis.text.y = element_text(size = 12),  # Adjust y-axis text size
    axis.title.x = element_text(size = 14),  # Adjust x-axis title size
    axis.title.y = element_text(size = 14),  # Adjust y-axis title size
    plot.title = element_text(size = 16),  # Adjust plot title size
    legend.text = element_text(size = 10),  # Adjust legend text size
    legend.title = element_text(size = 12)  # Adjust legend title size
    ) +

# Place the legend for genus colors on the right
theme(legend.position = "bottom")

ggsave(".tiff")
ggsave("16x12.pdf", width = 16, height = 12)



























