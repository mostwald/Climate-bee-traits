# Bee Pilosity and Lightness Analysis
# February 12, 2024

#### Set-up ####

#Load necessary packages
library(car) 
library(sf)
library(ggplot2)
library(maps)
library(lme4)
library(raster)
library(ape)
library(phytools)
library(geiger)
library(caper)
library(tidyverse)
library(nlme)
library(dplyr)
library(MuMIn) # to get R2 from GLMM
library(ggtree)
library(ggnewscale)
library(viridis)
library(emmeans)
library(scales)
library(paletteer)
library(RColorBrewer)



# Load in dataframe
setwd("~/Dropbox/UCSB/Manuscripts/Pilosity ML")
df <- read.csv("merged_df_jun12.csv") 

# Organize data
df$abs_latitude <- abs(df$Latitude)
df <- df[complete.cases(df$Longitude, df$Latitude), ]# filter for only specimens with geographic info
df <- df[complete.cases(df$pilosity_percentage), ] # filter for only specimens with pilosity data
df <- df[complete.cases(df$Brightness.Median), ] # filter for only specimens with lightness data

#### Descriptive statistics ####
aggregate(x = df$pilosity_percentage,               
          by = list(df$Family),             
          FUN = function(x) mean(x, na.rm = TRUE))  
mean(df$pilosity_percentage, na.rm=TRUE)
sd(df$pilosity_percentage, na.rm=TRUE)
range(df$pilosity_percentage, na.rm=TRUE)


aggregate(x = df$Brightness.Median,             
          by = list(df$Family),            
          FUN = function(x) mean(x, na.rm = TRUE))  
mean(df$Brightness.Median, na.rm=TRUE)
sd(df$Brightness.Median, na.rm=TRUE)
range(df$Brightness.Median, na.rm=TRUE)

# Fit mixed-effects models to estimate variance components
pilosity_variance <- lmer(pilosity_percentage ~ (1|Family) + (1|Subfamily) + (1|Genus), data = df)
pilosity_variance <- lmer(pilosity_percentage ~ (1|Family/Subfamily/Genus), data = df)
summary(pilosity_variance)
lightness_variance <- lmer(Brightness.Median ~ (1|Family) + (1|Subfamily) + (1|Genus), data = df)
summary(lightness_variance)

# Extract variance components
pilosity_var_comp <- as.data.frame(VarCorr(pilosity_variance))
print(pilosity_var_comp)
lightness_var_comp <- as.data.frame(VarCorr(lightness_variance))
print(lightness_var_comp)


#### Correlation test between pilosity and lightness ####
correlation_model <- lm(pilosity_percentage ~ Brightness.Median, data = df)
# Check assumptions
qqPlot(residuals(correlation_model)) 
plot(residuals(correlation_model) ~ fitted(correlation_model)) 

cor.test(df$pilosity_percentage, df$Brightness.Median,  method = "pearson")
#cor.test(df$pilosity_percentage, df$Brightness.Median,  method = "spearman")



#### Read in phylogenetic data ####

tree= read.tree(file="BEE_mat7gen_p8pmAa_fst.nwk") # read in tree

# Add missing taxa to most closely related genera
closest_genus_Belliturgula <- which(tree$tip.label == "Flavomeliturgula_tapana~Andrenidae~Panurginae~Panurgini")
tree_tips_added <- bind.tip(tree, "Belliturgula", where = closest_genus_Belliturgula)
tree_tips_added$edge.length[is.na(tree_tips_added$edge.length) | tree_tips_added$edge.length <= 1e-12] <- 0.01 # Ensure branch lengths are not NA or too small
closest_genus_Khuzimelissa <- which(tree$tip.label == "Flavomeliturgula_tapana~Andrenidae~Panurginae~Panurgini")
tree_tips_added <- bind.tip(tree_tips_added, "Khuzimelissa", where = closest_genus_Khuzimelissa)
tree_tips_added$edge.length[is.na(tree_tips_added$edge.length) | tree_tips_added$edge.length <= 1e-12] <- 0.01
closest_genus_Litocalliopsis <- which(tree$tip.label == "Calliopsis_anthidia~Andrenidae~Panurginae~Calliopsini") 
tree_tips_added <- bind.tip(tree_tips_added, "Litocalliopsis", where = closest_genus_Litocalliopsis)
tree_tips_added$edge.length[is.na(tree_tips_added$edge.length) | tree_tips_added$edge.length <= 1e-12] <- 0.01
closest_genus_Pseudopanurgus <- which(tree$tip.label == "Protandrena_verbesinae~Andrenidae~Panurginae~Protandrenini")
tree_tips_added <- bind.tip(tree_tips_added, "Pseudopanurgus", where = closest_genus_Pseudopanurgus)
tree_tips_added$edge.length[is.na(tree_tips_added$edge.length) | tree_tips_added$edge.length <= 1e-12] <- 0.01
closest_genus_Aparatrigona <- which(tree$tip.label == "Paratrigona_pacifica~Apidae~Apinae~Meliponini")
tree_tips_added <- bind.tip(tree_tips_added, "Aparatrigona", where = closest_genus_Aparatrigona)
tree_tips_added$edge.length[is.na(tree_tips_added$edge.length) | tree_tips_added$edge.length <= 1e-12] <- 0.01
closest_genus_Camargoia <- which(tree$tip.label == "Ptilotrigona_pereneae~Apidae~Apinae~Meliponini")
tree_tips_added <- bind.tip(tree_tips_added, "Camargoia", where = closest_genus_Camargoia)
tree_tips_added$edge.length[is.na(tree_tips_added$edge.length) | tree_tips_added$edge.length <= 1e-12] <- 0.01
closest_genus_Cemolobus <- which(tree$tip.label == "Eucera_nigrilabris~Apidae~Eucerinae~Eucerini")
tree_tips_added <- bind.tip(tree_tips_added, "Cemolobus", where = closest_genus_Cemolobus)
tree_tips_added$edge.length[is.na(tree_tips_added$edge.length) | tree_tips_added$edge.length <= 1e-12] <- 0.01
closest_genus_Cleptotrigona <- which(tree$tip.label == "Liotrigona_madecassa~Apidae~Apinae~Meliponini")
tree_tips_added <- bind.tip(tree_tips_added, "Cleptotrigona", where = closest_genus_Cleptotrigona)
tree_tips_added$edge.length[is.na(tree_tips_added$edge.length) | tree_tips_added$edge.length <= 1e-12] <- 0.01
closest_genus_Deltoptilla <- which(tree$tip.label == "Habropoda_laboriosa~Apidae~Anthophorinae~Anthophorini")
tree_tips_added <- bind.tip(tree_tips_added, "Deltoptilla", where = closest_genus_Deltoptilla)
tree_tips_added$edge.length[is.na(tree_tips_added$edge.length) | tree_tips_added$edge.length <= 1e-12] <- 0.01
closest_genus_Lophotrigona <- which(tree$tip.label == "Homotrigona_haematoptera~Apidae~Apinae~Meliponini")
tree_tips_added <- bind.tip(tree_tips_added, "Lophotrigona", where = closest_genus_Lophotrigona)
tree_tips_added$edge.length[is.na(tree_tips_added$edge.length) | tree_tips_added$edge.length <= 1e-12] <- 0.01
closest_genus_Papuatrigona <- which(tree$tip.label == "Lepidotrigona_terminata~Apidae~Apinae~Meliponini")
tree_tips_added <- bind.tip(tree_tips_added, "Papuatrigona", where = closest_genus_Papuatrigona)
tree_tips_added$edge.length[is.na(tree_tips_added$edge.length) | tree_tips_added$edge.length <= 1e-12] <- 0.01
closest_genus_Parapartomona <- which(tree$tip.label == "Partamona_testacea~Apidae~Apinae~Meliponini")
tree_tips_added <- bind.tip(tree_tips_added, "Parapartomona", where = closest_genus_Parapartomona)
tree_tips_added$edge.length[is.na(tree_tips_added$edge.length) | tree_tips_added$edge.length <= 1e-12] <- 0.01
closest_genus_Peponapis <- which(tree$tip.label == "Eucera_nigrilabris~Apidae~Eucerinae~Eucerini")
tree_tips_added <- bind.tip(tree_tips_added, "Peponapis", where = closest_genus_Peponapis)
tree_tips_added$edge.length[is.na(tree_tips_added$edge.length) | tree_tips_added$edge.length <= 1e-12] <- 0.01
closest_genus_Syntrichalonia <- which(tree$tip.label == "Eucera_nigrilabris~Apidae~Eucerinae~Eucerini")
tree_tips_added <- bind.tip(tree_tips_added, "Syntrichalonia", where = closest_genus_Syntrichalonia)
tree_tips_added$edge.length[is.na(tree_tips_added$edge.length) | tree_tips_added$edge.length <= 1e-12] <- 0.01
closest_genus_Tetraloniella <- which(tree$tip.label == "Eucera_nigrilabris~Apidae~Eucerinae~Eucerini")
tree_tips_added <- bind.tip(tree_tips_added, "Tetraloniella", where = closest_genus_Tetraloniella)
tree_tips_added$edge.length[is.na(tree_tips_added$edge.length) | tree_tips_added$edge.length <= 1e-12] <- 0.01
closest_genus_Pseudepeolus <- which(tree$tip.label == "Epeolus_scutellaris~Apidae~Nomadinae~Epeolini")
tree_tips_added <- bind.tip(tree_tips_added, "Pseudepeolus", where = closest_genus_Pseudepeolus)
tree_tips_added$edge.length[is.na(tree_tips_added$edge.length) | tree_tips_added$edge.length <= 1e-12] <- 0.01
closest_genus_Crawfordapis <- which(tree$tip.label == "Euryglossina_globuliceps~Colletidae~Euryglossinae~na")
tree_tips_added <- bind.tip(tree_tips_added, "Crawfordapis", where = closest_genus_Crawfordapis)
tree_tips_added$edge.length[is.na(tree_tips_added$edge.length) | tree_tips_added$edge.length <= 1e-12] <- 0.01
closest_genus_Hemirhiza <- which(tree$tip.label == "Meroglossa_itamuca~Colletidae~Hylaeinae~na")
tree_tips_added <- bind.tip(tree_tips_added, "Hemirhiza", where = closest_genus_Hemirhiza)
tree_tips_added$edge.length[is.na(tree_tips_added$edge.length) | tree_tips_added$edge.length <= 1e-12] <- 0.01
closest_genus_Agapostemonoides <- which(tree$tip.label == "Agapostemon_tyleri~Halictidae~Halictinae~Caenohalictini")
tree_tips_added <- bind.tip(tree_tips_added, "Agapostemonoides", where = closest_genus_Agapostemonoides)
tree_tips_added$edge.length[is.na(tree_tips_added$edge.length) | tree_tips_added$edge.length <= 1e-12] <- 0.01
closest_genus_Chlerogella <- which(tree$tip.label == "Ischnomelissa_zonata~Halictidae~Halictinae~Augochlorini")
tree_tips_added <- bind.tip(tree_tips_added, "Chlerogella", where = closest_genus_Chlerogella)
tree_tips_added$edge.length[is.na(tree_tips_added$edge.length) | tree_tips_added$edge.length <= 1e-12] <- 0.01
closest_genus_Chlerogelloides <- which(tree$tip.label == "Ischnomelissa_zonata~Halictidae~Halictinae~Augochlorini")
tree_tips_added <- bind.tip(tree_tips_added, "Chlerogelloides", where = closest_genus_Chlerogelloides)
tree_tips_added$edge.length[is.na(tree_tips_added$edge.length) | tree_tips_added$edge.length <= 1e-12] <- 0.01
closest_genus_Dinagapostemon <- which(tree$tip.label == "Agapostemon_tyleri~Halictidae~Halictinae~Caenohalictini")
tree_tips_added <- bind.tip(tree_tips_added, "Dinagapostemon", where = closest_genus_Dinagapostemon)
tree_tips_added$edge.length[is.na(tree_tips_added$edge.length) | tree_tips_added$edge.length <= 1e-12] <- 0.01
closest_genus_Paragapostemon <- which(tree$tip.label == "Agapostemon_tyleri~Halictidae~Halictinae~Caenohalictini")
tree_tips_added <- bind.tip(tree_tips_added, "Paragapostemon", where = closest_genus_Paragapostemon)
tree_tips_added$edge.length[is.na(tree_tips_added$edge.length) | tree_tips_added$edge.length <= 1e-12] <- 0.01
closest_genus_Pseudoagapostemon <- which(tree$tip.label == "Agapostemon_tyleri~Halictidae~Halictinae~Caenohalictini")
tree_tips_added <- bind.tip(tree_tips_added, "Pseudoagapostemon", where = closest_genus_Pseudoagapostemon)
tree_tips_added$edge.length[is.na(tree_tips_added$edge.length) | tree_tips_added$edge.length <= 1e-12] <- 0.01
closest_genus_Euptersia <- which(tree$tip.label == "Sphecodes_autumnalis~Halictidae~Halictinae~Sphecodini")
tree_tips_added <- bind.tip(tree_tips_added, "Euptersia", where = closest_genus_Euptersia)
tree_tips_added$edge.length[is.na(tree_tips_added$edge.length) | tree_tips_added$edge.length <= 1e-12] <- 0.01
closest_genus_Microsphecodes <- which(tree$tip.label == "Sphecodes_autumnalis~Halictidae~Halictinae~Sphecodini")
tree_tips_added <- bind.tip(tree_tips_added, "Microsphecodes", where = closest_genus_Microsphecodes)
tree_tips_added$edge.length[is.na(tree_tips_added$edge.length) | tree_tips_added$edge.length <= 1e-12] <- 0.01
closest_genus_Nesosphecodes <- which(tree$tip.label == "Sphecodes_autumnalis~Halictidae~Halictinae~Sphecodini")
tree_tips_added <- bind.tip(tree_tips_added, "Nesosphecodes", where = closest_genus_Nesosphecodes)
tree_tips_added$edge.length[is.na(tree_tips_added$edge.length) | tree_tips_added$edge.length <= 1e-12] <- 0.01
closest_genus_Habralictus <- which(tree$tip.label == "Caenohalictus_mourei~Halictidae~Halictinae~Caenohalictini")
tree_tips_added <- bind.tip(tree_tips_added, "Habralictus", where = closest_genus_Habralictus)
tree_tips_added$edge.length[is.na(tree_tips_added$edge.length) | tree_tips_added$edge.length <= 1e-12] <- 0.01
closest_genus_Halictonomia <- which(tree$tip.label == "Nomia_melanderi~Halictidae~Nomiinae~na")
tree_tips_added <- bind.tip(tree_tips_added, "Halictonomia", where = closest_genus_Halictonomia)
tree_tips_added$edge.length[is.na(tree_tips_added$edge.length) | tree_tips_added$edge.length <= 1e-12] <- 0.01
closest_genus_Melittidia <- which(tree$tip.label == "Nomia_melanderi~Halictidae~Nomiinae~na")
tree_tips_added <- bind.tip(tree_tips_added, "Melittidia", where = closest_genus_Melittidia)
tree_tips_added$edge.length[is.na(tree_tips_added$edge.length) | tree_tips_added$edge.length <= 1e-12] <- 0.01
closest_genus_Ptilonomia <- which(tree$tip.label == "Nomia_melanderi~Halictidae~Nomiinae~na")
tree_tips_added <- bind.tip(tree_tips_added, "Ptilonomia", where = closest_genus_Ptilonomia)
tree_tips_added$edge.length[is.na(tree_tips_added$edge.length) | tree_tips_added$edge.length <= 1e-12] <- 0.01
closest_genus_Reepenia <- which(tree$tip.label == "Nomia_melanderi~Halictidae~Nomiinae~na")
tree_tips_added <- bind.tip(tree_tips_added, "Reepenia", where = closest_genus_Reepenia)
tree_tips_added$edge.length[is.na(tree_tips_added$edge.length) | tree_tips_added$edge.length <= 1e-12] <- 0.01
closest_genus_Sphegocephala <- which(tree$tip.label == "Nomia_melanderi~Halictidae~Nomiinae~na")
tree_tips_added <- bind.tip(tree_tips_added, "Sphegocephala", where = closest_genus_Sphegocephala)
tree_tips_added$edge.length[is.na(tree_tips_added$edge.length) | tree_tips_added$edge.length <= 1e-12] <- 0.01
closest_genus_Megaloptina <- which(tree$tip.label == "Megommation_insigne~Halictidae~Halictinae~Augochlorini")
tree_tips_added <- bind.tip(tree_tips_added, "Megaloptina", where = closest_genus_Megaloptina)
tree_tips_added$edge.length[is.na(tree_tips_added$edge.length) | tree_tips_added$edge.length <= 1e-12] <- 0.01
closest_genus_Trichommation <- which(tree$tip.label == "Megommation_insigne~Halictidae~Halictinae~Augochlorini")
tree_tips_added <- bind.tip(tree_tips_added, "Trichommation", where = closest_genus_Trichommation)
tree_tips_added$edge.length[is.na(tree_tips_added$edge.length) | tree_tips_added$edge.length <= 1e-12] <- 0.01
closest_genus_Micralictoides <- which(tree$tip.label == "Dufourea_novaeangliae~Halictidae~Rophitinae~Rophitini")
tree_tips_added <- bind.tip(tree_tips_added, "Micralictoides", where = closest_genus_Micralictoides)
tree_tips_added$edge.length[is.na(tree_tips_added$edge.length) | tree_tips_added$edge.length <= 1e-12] <- 0.01
closest_genus_Morawitzia <- which(tree$tip.label == "Rophites_algirus~Halictidae~Rophitinae~Rophitini")
tree_tips_added <- bind.tip(tree_tips_added, "Morawitzia", where = closest_genus_Morawitzia)
tree_tips_added$edge.length[is.na(tree_tips_added$edge.length) | tree_tips_added$edge.length <= 1e-12] <- 0.01
closest_genus_Afrostelis <- which(tree$tip.label == "Stelis_punctulatissima~Megachilidae~Megachilinae~Anthidiini")
tree_tips_added <- bind.tip(tree_tips_added, "Afrostelis", where = closest_genus_Afrostelis)
tree_tips_added$edge.length[is.na(tree_tips_added$edge.length) | tree_tips_added$edge.length <= 1e-12] <- 0.01
closest_genus_Apianthidium <- which(tree$tip.label == "Trachusa_larreae~Megachilidae~Megachilinae~Anthidiini")
tree_tips_added <- bind.tip(tree_tips_added, "Apianthidium", where = closest_genus_Apianthidium)
tree_tips_added$edge.length[is.na(tree_tips_added$edge.length) | tree_tips_added$edge.length <= 1e-12] <- 0.01
closest_genus_Hypanthidioides <- which(tree$tip.label == "Anthodioctes_mapirensis~Megachilidae~Megachilinae~Anthidiini")
tree_tips_added <- bind.tip(tree_tips_added, "Hypanthidioides", where = closest_genus_Hypanthidioides)
tree_tips_added$edge.length[is.na(tree_tips_added$edge.length) | tree_tips_added$edge.length <= 1e-12] <- 0.01
closest_genus_Microthurge <- which(tree$tip.label == "Lithurgus_pullatus~Megachilidae~Lithurginae~Lithurgini")
tree_tips_added <- bind.tip(tree_tips_added, "Microthurge", where = closest_genus_Microthurge)
tree_tips_added$edge.length[is.na(tree_tips_added$edge.length) | tree_tips_added$edge.length <= 1e-12] <- 0.01
closest_genus_Noteriades <- which(tree$tip.label == "Megachile_rotundata~Megachilidae~Megachilinae~Megachilini")
tree_tips_added <- bind.tip(tree_tips_added, "Noteriades", where = closest_genus_Noteriades)
tree_tips_added$edge.length[is.na(tree_tips_added$edge.length) | tree_tips_added$edge.length <= 1e-12] <- 0.01
closest_genus_Xenostelis <- which(tree$tip.label == "Plesianthidium_rufocaudatum~Megachilidae~Megachilinae~Anthidiini")
tree_tips_added <- bind.tip(tree_tips_added, "Xenostelis", where = closest_genus_Xenostelis)
tree_tips_added$edge.length[is.na(tree_tips_added$edge.length) | tree_tips_added$edge.length <= 1e-12] <- 0.01
closest_genus_Ceratomonia <- which(tree$tip.label == "Meganomia_binghami~Melittidae~Meganomiinae_na~")
tree_tips_added <- bind.tip(tree_tips_added, "Ceratomonia", where = closest_genus_Ceratomonia)
tree_tips_added$edge.length[is.na(tree_tips_added$edge.length) | tree_tips_added$edge.length <= 1e-12] <- 0.01
closest_genus_Pseudophilanthus <- which(tree$tip.label == "Meganomia_binghami~Melittidae~Meganomiinae_na~")
tree_tips_added <- bind.tip(tree_tips_added, "Pseudophilanthus", where = closest_genus_Pseudophilanthus)
tree_tips_added$edge.length[is.na(tree_tips_added$edge.length) | tree_tips_added$edge.length <= 1e-12] <- 0.01
closest_genus_Alepidoscelis <- which(tree$tip.label == "Ptilothrix_bombiformis~Apidae~Eucerinae~Emphorini")
tree_tips_added <- bind.tip(tree_tips_added, "Alepidoscelis", where = closest_genus_Alepidoscelis)
tree_tips_added$edge.length[is.na(tree_tips_added$edge.length) | tree_tips_added$edge.length <= 1e-12] <- 0.01
closest_genus_Xeranthrena <- which(tree$tip.label == "Callonychium_atacamense~Andrenidae~Panurginae~Calliopsini")
tree_tips_added <- bind.tip(tree_tips_added, "Xeranthrena", where = closest_genus_Xeranthrena)
tree_tips_added$edge.length[is.na(tree_tips_added$edge.length) | tree_tips_added$edge.length <= 1e-12] <- 0.01

# Extract genus from tree tip labels
genus_from_tree_tips_added <- sub("_.*", "", tree_tips_added$tip.label)

# Match genus with df$Genus
matching_indices_tips_added <- match(df$Genus, genus_from_tree_tips_added)

# Create a new column Taxonomic_Info in df
df$Taxonomic_Info <- tree_tips_added$tip.label[matching_indices_tips_added]

# Extract tip labels from Taxonomic_Info column
tip_labels_to_keep_tips_added <- df$Taxonomic_Info

# Prune the tree based on the extracted tip labels
pruned_tree <- drop.tip(tree_tips_added, setdiff(tree_tips_added$tip.label, tip_labels_to_keep_tips_added))

# Extract tip labels from pruned_tree
pruned_tip_labels<- pruned_tree$tip.label

# Filter rows in df that have corresponding tips in pruned_tree
matched_df <- df[df$Taxonomic_Info %in% pruned_tip_labels, ]

## pruned_tree_tips_added has an additional 60 genera (319 vs 369)






#### Is there phylogenetic signal in pilosity and lightness? ####

## Aggregate data at the genus level
mean_pilosity_phylosig <- aggregate(pilosity_percentage ~ Taxonomic_Info, data = matched_df, FUN = mean)
mean_pilosity_vector <- setNames(mean_pilosity_phylosig$pilosity_percentage, mean_pilosity_phylosig$Taxonomic_Info)
phylosig(pruned_tree, mean_pilosity_vector,method="lambda",test=TRUE)

mean_lightness_phylosig <- aggregate(Brightness.Median ~ Taxonomic_Info, data = matched_df, FUN = mean)
mean_lightness_vector <- setNames(mean_lightness_phylosig$Brightness.Median, mean_lightness_phylosig$Taxonomic_Info)
phylosig(pruned_tree, mean_lightness_vector,method="lambda",test=TRUE)







#### Bioclimatic Data ####

# Import WorldClim Mean Monthly Precip data

prec1 <- raster("wc2.1_30s_prec_01.tif")
prec2 <- raster("wc2.1_30s_prec_02.tif")
prec3 <- raster("wc2.1_30s_prec_03.tif")
prec4 <- raster("wc2.1_30s_prec_04.tif")
prec5 <- raster("wc2.1_30s_prec_05.tif")
prec6 <- raster("wc2.1_30s_prec_06.tif")
prec7 <- raster("wc2.1_30s_prec_07.tif")
prec8 <- raster("wc2.1_30s_prec_08.tif")
prec9 <- raster("wc2.1_30s_prec_09.tif")
prec10 <- raster("wc2.1_30s_prec_10.tif")
prec11 <- raster("wc2.1_30s_prec_11.tif")
prec12 <- raster("wc2.1_30s_prec_12.tif")

# Import WorldClim Mean Monthly Vapor pressure data

vapr1 <- raster("wc2.1_30s_vapr_01.tif")
vapr2 <- raster("wc2.1_30s_vapr_02.tif")
vapr3 <- raster("wc2.1_30s_vapr_03.tif")
vapr4 <- raster("wc2.1_30s_vapr_04.tif")
vapr5 <- raster("wc2.1_30s_vapr_05.tif")
vapr6 <- raster("wc2.1_30s_vapr_06.tif")
vapr7 <- raster("wc2.1_30s_vapr_07.tif")
vapr8 <- raster("wc2.1_30s_vapr_08.tif")
vapr9 <- raster("wc2.1_30s_vapr_09.tif")
vapr10 <- raster("wc2.1_30s_vapr_10.tif")
vapr11 <- raster("wc2.1_30s_vapr_11.tif")
vapr12 <- raster("wc2.1_30s_vapr_12.tif")

# Import WorldClim Mean Monthly Solar Radiation data

srad1 <- raster("wc2.1_30s_srad_01.tif")
srad2 <- raster("wc2.1_30s_srad_02.tif")
srad3 <- raster("wc2.1_30s_srad_03.tif")
srad4 <- raster("wc2.1_30s_srad_04.tif")
srad5 <- raster("wc2.1_30s_srad_05.tif")
srad6 <- raster("wc2.1_30s_srad_06.tif")
srad7 <- raster("wc2.1_30s_srad_07.tif")
srad8 <- raster("wc2.1_30s_srad_08.tif")
srad9 <- raster("wc2.1_30s_srad_09.tif")
srad10 <- raster("wc2.1_30s_srad_10.tif")
srad11 <- raster("wc2.1_30s_srad_11.tif")
srad12 <- raster("wc2.1_30s_srad_12.tif")

# Import WorldClim Mean Monthly Min Temp data

tmin1 <- raster("wc2.1_30s_tmin_01.tif")
tmin2 <- raster("wc2.1_30s_tmin_02.tif")
tmin3 <- raster("wc2.1_30s_tmin_03.tif")
tmin4 <- raster("wc2.1_30s_tmin_04.tif")
tmin5 <- raster("wc2.1_30s_tmin_05.tif")
tmin6 <- raster("wc2.1_30s_tmin_06.tif")
tmin7 <- raster("wc2.1_30s_tmin_07.tif")
tmin8 <- raster("wc2.1_30s_tmin_08.tif")
tmin9 <- raster("wc2.1_30s_tmin_09.tif")
tmin10 <- raster("wc2.1_30s_tmin_10.tif")
tmin11 <- raster("wc2.1_30s_tmin_11.tif")
tmin12 <- raster("wc2.1_30s_tmin_12.tif")

# Import WorldClim Mean Monthly Avg Temp data

tavg1 <- raster("wc2.1_30s_tavg_01.tif")
tavg2 <- raster("wc2.1_30s_tavg_02.tif")
tavg3 <- raster("wc2.1_30s_tavg_03.tif")
tavg4 <- raster("wc2.1_30s_tavg_04.tif")
tavg5 <- raster("wc2.1_30s_tavg_05.tif")
tavg6 <- raster("wc2.1_30s_tavg_06.tif")
tavg7 <- raster("wc2.1_30s_tavg_07.tif")
tavg8 <- raster("wc2.1_30s_tavg_08.tif")
tavg9 <- raster("wc2.1_30s_tavg_09.tif")
tavg10 <- raster("wc2.1_30s_tavg_10.tif")
tavg11 <- raster("wc2.1_30s_tavg_11.tif")
tavg12 <- raster("wc2.1_30s_tavg_12.tif")

#Import WorldClim Mean Monthly max Temp data

tmax1 <- raster("wc2.1_30s_tmax_01.tif")
tmax2 <- raster("wc2.1_30s_tmax_02.tif")
tmax3 <- raster("wc2.1_30s_tmax_03.tif")
tmax4 <- raster("wc2.1_30s_tmax_04.tif")
tmax5 <- raster("wc2.1_30s_tmax_05.tif")
tmax6 <- raster("wc2.1_30s_tmax_06.tif")
tmax7 <- raster("wc2.1_30s_tmax_07.tif")
tmax8 <- raster("wc2.1_30s_tmax_08.tif")
tmax9 <- raster("wc2.1_30s_tmax_09.tif")
tmax10 <- raster("wc2.1_30s_tmax_10.tif")
tmax11 <- raster("wc2.1_30s_tmax_11.tif")
tmax12 <- raster("wc2.1_30s_tmax_12.tif")


# Create lists to store the climate rasters
precip_rasters <- list(prec1, prec2, prec3, prec4, prec5, prec6, prec7, prec8, prec9, prec10, prec11, prec12)
vapr_rasters <- list(vapr1, vapr2, vapr3, vapr4, vapr5, vapr6, vapr7, vapr8, vapr9, vapr10, vapr11, vapr12)
srad_rasters <- list(srad1, srad2, srad3, srad4, srad5, srad6, srad7, srad8, srad9, srad10, srad11, srad12)
tmin_rasters <- list(tmin1, tmin2, tmin3, tmin4, tmin5, tmin6, tmin7, tmin8, tmin9, tmin10, tmin11, tmin12)
tavg_rasters <- list(tavg1, tavg2, tavg3, tavg4, tavg5, tavg6, tavg7, tavg8, tavg9, tavg10, tavg11, tavg12)
tmax_rasters <- list(tmax1, tmax2, tmax3, tmax4, tmax5, tmax6, tmax7, tmax8, tmax9, tmax10, tmax11, tmax12)

# Create empty matrices to store climate values for each location
precip_values <- matrix(NA, nrow = nrow(df), ncol = 12) 
vapr_values <- matrix(NA, nrow = nrow(df), ncol = 12)  
srad_values <- matrix(NA, nrow = nrow(df), ncol = 12)  
tmin_values <- matrix(NA, nrow = nrow(df), ncol = 12)
tavg_values <- matrix(NA, nrow = nrow(df), ncol = 12) 
tmax_values <- matrix(NA, nrow = nrow(df), ncol = 12)  

# Loop through each month to extract climate values for each location
for (i in 1:12) {precip_values[, i] <- raster::extract(precip_rasters[[i]], df[, c("Longitude", "Latitude")]) }
for (i in 1:12) {vapr_values[, i] <- raster::extract(vapr_rasters[[i]], df[, c("Longitude", "Latitude")])}
for (i in 1:12) {srad_values[, i] <- raster::extract(srad_rasters[[i]], df[, c("Longitude", "Latitude")])}
for (i in 1:12) {tmin_values[, i] <- raster::extract(tmin_rasters[[i]], df[, c("Longitude", "Latitude")])} 
for (i in 1:12) {tavg_values[, i] <- raster::extract(tavg_rasters[[i]], df[, c("Longitude", "Latitude")])}
for (i in 1:12) {tmax_values[, i] <- raster::extract(tmax_rasters[[i]], df[, c("Longitude", "Latitude")])}

# Calculate climate variable values
df$precip_annual <- rowSums(precip_values, na.rm = TRUE)
df$precip_seasonality <- (apply(precip_values, 1, sd, na.rm = TRUE))*10
df$vapr <- rowMeans(vapr_values, na.rm = TRUE)
df$srad <- (rowMeans(srad_values, na.rm = TRUE))/1000
df$tmin_coldest_month <- apply(tmin_values, 1, max, na.rm = TRUE)
df$temp_annual <- rowMeans(tavg_values, na.rm = TRUE)
df$temp_seasonality <- (apply(tavg_values, 1, sd, na.rm = TRUE))*100
df$tmax_hottest_month <- apply(tmax_values, 1, max, na.rm = TRUE)

df <- df[complete.cases(df$srad), ] # filter for only specimens with complete bioclimatic data
view(df)
  
  

#### Assess normality of predictors and response variables ####

# Add transformations to dataset
df$sqrt_precip_annual <- sqrt(df$precip_annual)
df$log_precip_seasonality <- log(df$precip_seasonality)
df$sqrt_temp_seasonality <- sqrt(df$temp_seasonality)
df$sqrt_lightness <- sqrt(df$Brightness.Median)


hist(df$precip_annual)
hist(df$sqrt_precip_annual) # better
hist(df$precip_seasonality) 
hist(df$log_precip_seasonality) # better
hist(df$srad) 
hist(df$vapr) 
hist(df$temp_annual) 
hist(df$temp_seasonality)
hist(sqrt(df$temp_seasonality)) # better
hist(df$tmin_coldest_month) 
hist(df$tmax_hottest_month) 
hist(df$Brightness.Mean)
hist(df$pilosity_percentage)
hist(df$Brightness.Median)
hist(sqrt(df$Brightness.Median)) # better

qqPlot(df$precip_annual)
qqPlot(df$sqrt_precip_annual) # better
qqPlot(df$precip_seasonality) 
qqPlot(df$log_precip_seasonality) # better
qqPlot(df$srad) 
qqPlot(df$vapr) 
qqPlot(df$temp_annual) 
qqPlot(df$temp_seasonality)
qqPlot(sqrt(df$temp_seasonality)) # better
qqPlot(df$tmin_coldest_month) 
qqPlot(df$tmax_hottest_month) 
qqPlot(df$Brightness.Median)
qqPlot(sqrt(df$Brightness.Median))
qqPlot(df$pilosity_percentage)




#### Do bioclimatic variables predict trait variation? ##### 

# First, to make sure we're not biasing our analysis with the order of removal of variables, test each predictor separately (see Kellerman 2018)
pilosity_glmm_sqrt_precip_annual <- lmer(pilosity_percentage ~ sqrt_precip_annual+ (1|Subfamily/Genus), data = df)
Anova(pilosity_glmm_sqrt_precip_annual, type="III") # P = 0.02673 *  
r.squaredGLMM(pilosity_glmm_sqrt_precip_annual) # 0.008820317 # look at marginal R2 to look at fixed effects

pilosity_glmm_log_precip_seasonality <- lmer(pilosity_percentage ~ log_precip_seasonality+ (1|Subfamily/Genus), data = df)
Anova(pilosity_glmm_log_precip_seasonality, type="III") # P = 0.2183    

pilosity_glmm_temp_annual <- lmer(pilosity_percentage ~ temp_annual+ (1|Subfamily/Genus), data = df)
Anova(pilosity_glmm_temp_annual, type="III") # P = 0.01397 *  
r.squaredGLMM(pilosity_glmm_temp_annual) # 0.01086405

pilosity_glmm_sqrt_temp_seasonality <- lmer(pilosity_percentage ~ sqrt_temp_seasonality+ (1|Subfamily/Genus), data = df)
Anova(pilosity_glmm_sqrt_temp_seasonality, type="III") # P =  0.896  

pilosity_glmm_vapr <- lmer(pilosity_percentage ~ vapr+ (1|Subfamily/Genus), data = df)
Anova(pilosity_glmm_vapr, type="III") # P = 0.4815

pilosity_glmm_srad <- lmer(pilosity_percentage ~ srad + (1|Subfamily/Genus), data = df)
Anova(pilosity_glmm_srad, type="III") # P = 0.01719 *  
r.squaredGLMM(pilosity_glmm_srad) # 0.009050986

pilosity_glmm_tmax <- lmer(pilosity_percentage ~ tmax_hottest_month + (1|Subfamily/Genus), data = df)
Anova(pilosity_glmm_tmax, type="III") # P = 0.0003001 ***
r.squaredGLMM(pilosity_glmm_tmax) # 0.01967616

pilosity_glmm_tmin <- lmer(pilosity_percentage ~ tmin_coldest_month+ (1|Subfamily/Genus), data = df)
Anova(pilosity_glmm_tmin, type="III") # 0.01153 *  
r.squaredGLMM(pilosity_glmm_tmin) # 0.009347625

# Priority order of influential variables: tmax, temp_annual, tmin, srad, sqrt_precip_annual

pilosity_glmm1 <- lmer(pilosity_percentage ~ sqrt_precip_annual + log_precip_seasonality + temp_annual + sqrt_temp_seasonality + vapr + srad + tmin_coldest_month + tmax_hottest_month + (1|Subfamily/Genus), data = df)
vif(pilosity_glmm1) # remove sqrt_temp_seasonality 

pilosity_glmm2 <- lmer(pilosity_percentage ~ sqrt_precip_annual + log_precip_seasonality + tmin_coldest_month + temp_annual + vapr + srad + tmax_hottest_month + (1|Subfamily/Genus), data = df )
vif(pilosity_glmm2) # remove vapr

pilosity_glmm3 <- lmer(pilosity_percentage ~ sqrt_precip_annual + log_precip_seasonality + tmin_coldest_month + temp_annual + srad + tmax_hottest_month + (1|Subfamily/Genus), data =  df) 
vif(pilosity_glmm3) # remove tmin

pilosity_glmm4 <- lmer(pilosity_percentage ~ sqrt_precip_annual + log_precip_seasonality + temp_annual + srad + tmax_hottest_month + (1|Subfamily/Genus), data =  df) 
vif(pilosity_glmm4) # all VIF < 5
AIC(pilosity_glmm4) # 4601.81
Anova(pilosity_glmm4, type="III") # remove log_precip_seasonality

pilosity_glmm5 <- lmer(pilosity_percentage ~ sqrt_precip_annual + temp_annual + srad + tmax_hottest_month + (1|Subfamily/Genus), data =  df) 
vif(pilosity_glmm5) # all VIF < 5
AIC(pilosity_glmm5) # 4602.602
Anova(pilosity_glmm5, type="III") # remove srad

pilosity_glmm6 <- lmer(pilosity_percentage ~ sqrt_precip_annual + temp_annual + tmax_hottest_month + (1|Subfamily/Genus), data =  df) 
vif(pilosity_glmm6) # all VIF < 5
AIC(pilosity_glmm6) # 4600.664
Anova(pilosity_glmm6, type="III") # remove temp_annual
summary(pilosity_glmm6)

pilosity_glmm7 <- lmer(pilosity_percentage ~ sqrt_precip_annual + tmax_hottest_month + (1|Subfamily/Genus), data =  df) 
vif(pilosity_glmm7) # all VIF < 5
AIC(pilosity_glmm7) # 4599.417
Anova(pilosity_glmm7, type="III") # remove sqrt_precip_annual

pilosity_glmm8 <- lmer(pilosity_percentage ~ tmax_hottest_month + (1|Subfamily/Genus), data =  df) 
AIC(pilosity_glmm8) # 4596.064
summary(pilosity_glmm8)
r.squaredGLMM(pilosity_glmm8)
Anova(pilosity_glmm8, type="III")

# Check assumptions
qqPlot(residuals(pilosity_glmm8)) 
plot(residuals(pilosity_glmm8) ~ fitted(pilosity_glmm8))

# Check lm without random effect
pilosity_lm <- lm(pilosity_percentage ~ tmax_hottest_month, data = df)
AIC(pilosity_lm) # 4647.325

# Pilosity increases with tmax of the hottest month.

## Repeat this process for lightness (sqrt):

lightness_glmm_sqrt_precip_annual <- lmer(sqrt_lightness ~ sqrt_precip_annual+ (1|Subfamily/Genus), data = df)
Anova(lightness_glmm_sqrt_precip_annual, type="III") # 0.0007361 ***
r.squaredGLMM(lightness_glmm_sqrt_precip_annual) # 0.02152753 # look at marginal R2 to look at fixed effects

lightness_glmm_log_precip_seasonality <- lmer(sqrt_lightness ~ log_precip_seasonality+ (1|Subfamily/Genus), data = df)
Anova(lightness_glmm_log_precip_seasonality, type="III") # P = 0.03695 * 
r.squaredGLMM(lightness_glmm_log_precip_seasonality) # 0.008298748

lightness_glmm_temp_annual <- lmer(sqrt_lightness ~ temp_annual+ (1|Subfamily/Genus), data = df)
Anova(lightness_glmm_temp_annual, type="III") # P = 0.1433   

lightness_glmm_sqrt_temp_seasonality <- lmer(sqrt_lightness ~ sqrt_temp_seasonality+ (1|Subfamily/Genus), data = df)
Anova(lightness_glmm_sqrt_temp_seasonality, type="III") # 0.1088 

lightness_glmm_vapr <- lmer(sqrt_lightness ~ vapr+ (1|Subfamily/Genus), data = df)
Anova(lightness_glmm_vapr, type="III") # P = 0.5198   

lightness_glmm_srad <- lmer(sqrt_lightness ~ srad + (1|Subfamily/Genus), data = df)
Anova(lightness_glmm_srad, type="III") # P = 0.002699 ** 
r.squaredGLMM(lightness_glmm_srad) # 0.01527205

lightness_glmm_tmax <- lmer(sqrt_lightness ~ tmax_hottest_month + (1|Subfamily/Genus), data = df)
Anova(lightness_glmm_tmax, type="III") # P = 7.704e-05 ***
r.squaredGLMM(lightness_glmm_tmax) # 0.02494344

lightness_glmm_tmin <- lmer(sqrt_lightness ~ tmin_coldest_month+ (1|Subfamily/Genus), data = df)
Anova(lightness_glmm_tmin, type="III") # 0.02166 * 
r.squaredGLMM(lightness_glmm_tmin) # 0.008296559

# Priority order: tmax, sqrt_precip_annual, srad, tmin, log_precip_seasonality

lightness_glmm1 <- lmer(sqrt_lightness ~ sqrt_precip_annual + log_precip_seasonality + temp_annual + sqrt_temp_seasonality + vapr + srad + tmin_coldest_month + tmax_hottest_month + (1|Subfamily/Genus), data = df)
vif(lightness_glmm1) # remove temp_annual

lightness_glmm2 <- lmer(sqrt_lightness ~ sqrt_precip_annual + log_precip_seasonality + sqrt_temp_seasonality + vapr + srad + tmin_coldest_month + tmax_hottest_month + (1|Subfamily/Genus), data = df)
vif(lightness_glmm2) # remove vapr

lightness_glmm3 <- lmer(sqrt_lightness ~ sqrt_precip_annual + log_precip_seasonality + sqrt_temp_seasonality + srad + tmin_coldest_month + tmax_hottest_month + (1|Subfamily/Genus), data = df)
vif(lightness_glmm3) # all VIF < 5
AIC(lightness_glmm3) # 1900.366
Anova(lightness_glmm3, type="III") # remove temp_seasonality

lightness_glmm4 <- lmer(sqrt_lightness ~ sqrt_precip_annual + log_precip_seasonality + srad + tmin_coldest_month + tmax_hottest_month + (1|Subfamily/Genus), data = df)
AIC(lightness_glmm4) # 1891.522
Anova(lightness_glmm4, type="III") # remove precip_seasonality

lightness_glmm5 <- lmer(sqrt_lightness ~ sqrt_precip_annual + srad + tmin_coldest_month + tmax_hottest_month + (1|Subfamily/Genus), data = df)
AIC(lightness_glmm5) # 1886.349
Anova(lightness_glmm5, type="III") # remove tmin

lightness_glmm6 <- lmer(sqrt_lightness ~ sqrt_precip_annual + srad + tmax_hottest_month + (1|Subfamily/Genus), data = df)
AIC(lightness_glmm6) # 1879.934
Anova(lightness_glmm6, type="III") # remove srad

lightness_glmm7 <- lmer(sqrt_lightness ~ sqrt_precip_annual + tmax_hottest_month + (1|Subfamily/Genus), data = df)
AIC(lightness_glmm7) # 1874.609
Anova(lightness_glmm7, type="III")
r.squaredGLMM(lightness_glmm7)
summary(lightness_glmm7)

# lightness increases with tmax and decreases with sqrt_precip_annual.

# Check assumptions
qqPlot(residuals(lightness_glmm7)) 
plot(residuals(lightness_glmm7) ~ fitted(lightness_glmm7)) 


#### Does biome predict trait variation? ####

# Load in WWF ecoregion delims
wwf_biomes <- sf::st_read("wwf_terr_ecos.shp")

# Convert df to sf
sf <- st_as_sf(df, coords = c("Longitude", "Latitude"),crs = st_crs("+proj=longlat +datum=WGS84")) # Specify the CRS (coordinate reference system)

# check that coordinate projections are equal
st_crs(wwf_biomes) == st_crs(sf)

# Create datasets for each biome
biome_delim_occurrences <- list()

# Iterate over unique biome names
unique_biomes <- unique(wwf_biomes$BIOME)
for(i in 1:length(unique_biomes)){
  biome_shp <- wwf_biomes %>%
    filter(BIOME == unique_biomes[[i]]) %>%
    st_make_valid()
  
  sf_regioned <- sf[biome_shp, ]
  
  # Add biome ID to each of the subsets
  sf_regioned <- sf_regioned %>%
    mutate(wwf_biome = unique_biomes[[i]])
  
  # Append the subsetted data to the list
  biome_delim_occurrences[[i]] <- sf_regioned
}

# Create the output dataframe
df_biomes <- bind_rows(biome_delim_occurrences)
# Convert biome to a factor
df_biomes$wwf_biome <- factor(df_biomes$wwf_biome, levels = 1:14)
# Collapse levels
df_biomes <- df_biomes %>%
  mutate(collapsed_biomes = recode(wwf_biome,
                                   `1` = "Tropical & Subtropical Forests",
                                   `2` = "Tropical & Subtropical Forests",
                                   `3` = "Tropical & Subtropical Forests",
                                   `4` = "Temperate Forests",
                                   `5` = "Temperate Forests",
                                   `6` = "Misc.",
                                   `7` = "Grasslands & Shrublands",
                                   `8` = "Grasslands & Shrublands",
                                   `9` = "Grasslands & Shrublands",
                                   `10` = "Grasslands & Shrublands",
                                   `11` = "Misc.",
                                   `12` = "Temperate Forests",
                                   `13` = "Deserts & Xeric Shrublands",
                                   `14` = "Misc."))
# Biome key:
# 1: Tropical & Subtropical Moist Broadleaf Forests
# 2: Tropical & Subtropical Dry Broadleaf Forests
# 3: Tropical & Subtropical Coniferous Forests
# 4: Temperate Broadleaf & Mixed Forests
# 5: Temperate Conifer Forests
# 6: Boreal Forests/Taiga
# 7: Tropical & Subtropical Grasslands, Savannas & Shrublands
# 8: Temperate Grasslands, Savannas & Shrublands
# 9: Flooded Grasslands & Savannas
# 10: Montane Grasslands & Shrublands
# 11: Tundra
# 12: Mediterranean Forests, Woodlands & Scrub
# 13: Deserts & Xeric Shrublands
# 14: Mangroves

# Collapsed Biomes:
# 1, 2, 3: Tropical & Subtropical Forests
# 4, 5, 12: Temperate Forests
# 7, 8, 9, 10: Grasslands & Shrublands
# 13: Deserts & Xeric Shrublands
# 6, 11, 14: Misc

# Remove rows where biome is "Misc." or NA
biomes_subset <- subset(df_biomes, collapsed_biomes != "Misc." & !is.na(collapsed_biomes))

# Run GLMM on pilosity data
biome_glmm_pilosity <- lmer(pilosity_percentage ~ collapsed_biomes + (1|Subfamily/Genus), data =  biomes_subset)
summary(biome_glmm_pilosity)
AIC(biome_glmm_pilosity) # 4360.448
r.squaredGLMM(biome_glmm_pilosity)

# Check assumptions
qqPlot(residuals(biome_glmm_pilosity))
plot(residuals(biome_glmm_pilosity) ~ fitted(biome_glmm_pilosity)) 

# Compare to LM without random effect of taxonomic info
biome_pilosity_lm <- lm(pilosity_percentage ~ collapsed_biomes, data = biomes_subset)
AIC(biome_pilosity_lm) # 4413.971
r.squaredGLMM(biome_glmm_lightness)

# Run ANOVA and estimated marginal means
Anova(biome_glmm_pilosity, type="III")
emmeans_results_pilosity <- emmeans(biome_glmm_pilosity, ~ collapsed_biomes) # estimated marginal means for pairwise biome comparisons
pairwise_comparisons_pilosity <- pairs(emmeans_results_pilosity) # Perform pairwise comparisons between biomes
summary(pairwise_comparisons_pilosity)


# Now repeat for lightness
# Run GLMM on lightness data
biome_glmm_lightness <- lmer(Brightness.Median ~ collapsed_biomes + (1|Subfamily/Genus), data =  biomes_subset)
AIC(biome_glmm_lightness) # 4709.672

# Check assumptions
qqPlot(residuals(biome_glmm_lightness)) 
plot(residuals(biome_glmm_lightness) ~ fitted(biome_glmm_lightness)) 

# Compare to LM without random effect of taxonomic info
biome_lightness_lm <- lm(Brightness.Median ~ collapsed_biomes, data = biomes_subset)
AIC(biome_lightness_lm) # 4760.109

# Run ANOVA and estimated marginal means
Anova(biome_glmm_lightness, type="III")
emmeans_results_lightness <- emmeans(biome_glmm_lightness, ~ collapsed_biomes) # estimated marginal means for pairwise biome comparisons
pairwise_comparisons_lightness <- pairs(emmeans_results_lightness) # Perform pairwise comparisons between biomes
summary(pairwise_comparisons_lightness)





#### Plot biomes against trait values ####

library(cowplot)

pilosity_biome_plot <- ggplot(data = biomes_subset, aes(collapsed_biomes, pilosity_percentage)) +
  geom_boxplot(fill=c("#139D89", "#139D89", "#139D89","#F56B5C"),alpha=0.7)+
  labs(x = NULL,y = "Hair Coverage (%)", size=12) +
  theme_classic() +
  theme(axis.text.x = element_blank(),axis.text.y = element_text(color="black", size = 12), legend.position = "none") 

lightness_biome_plot <- ggplot(data = biomes_subset, aes(collapsed_biomes, Brightness.Median)) +
  geom_boxplot(fill=c("#139D89", "#139D89", "#139D89","#F56B5C"),alpha=0.7)+
  labs(x = NULL,y = "Lightness", size=12) +
  theme_classic() +
theme(axis.text.x = element_blank(),axis.text.y = element_text(color="black", size = 12), legend.position = "none")


#### Plot climate variables against trait values ####

family_palette <- c("#040313","#8707A6","#F56B5C","#FEBE2A","#139D89","#F9973E","#75D054")
# Plot pilosity vs. tmax
ggplot(data = df, aes(tmax_hottest_month, pilosity_percentage, color=Family)) +
  geom_point(alpha=0.8, size=3) +
  geom_smooth(method = "lm", color = "black") +
  labs(x = expression("Max. Temp. of the Hottest Month (" * degree * "C)"),y = "Hair Coverage (%)", size=12) +
  theme_classic() +
  scale_color_manual(values = family_palette) +
  theme(axis.text.x = element_text(color="black", size = 12),axis.text.y = element_text(color="black", size = 12))

# Plot lightness vs. tmax
ggplot(data = df, aes(tmax_hottest_month, sqrt_lightness, color=Family)) +
  geom_point(alpha=0.8, size=3) +
  geom_smooth(method = "lm", color = "black") +
  scale_color_manual(values = family_palette) +
  labs(x = expression("Max. Temp. of the Hottest Month (" * degree * "C)"),y = "Sqrt. Lightness", size=12) +
  theme_classic() +
  theme(axis.text.x = element_text(color="black", size = 12),axis.text.y = element_text(color="black", size = 12))


# Plot lightness vs. precip_annual
ggplot(data = df, aes(sqrt_precip_annual, sqrt_lightness, color=Family)) +
  geom_point(alpha=0.8, size=3) +
  geom_smooth(method = "lm", color = "black") +
  scale_color_manual(values = family_palette) +
  labs(x = expression("Sqrt. Annual Precipitation (mm)"),y = "Sqrt. Lightness", size=12) +
  theme_classic() +
  theme(axis.text.x = element_text(color="black", size = 12),axis.text.y = element_text(color="black", size = 12)) 


#### Plot pilosity vs. lightness ####

ggplot(data=df, aes(pilosity_percentage, Brightness.Median, color=Family)) +
  geom_point(alpha=0.5, size=2) +
  labs(y = "Lightness", x = "Hair Coverage (%)", size=12) + 
  scale_color_manual(values = family_palette) +
  theme_classic() +
  theme(axis.text.x = element_text(color="black", size = 12),axis.text.y = element_text(color="black", size = 12))


#### Plot world maps ####
world <- map_data("world")
ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group), fill = "gray", color = "gray", size = 0.1) +
  geom_point(data = df, aes(x = Longitude, y = Latitude, color = pilosity_percentage), size = 1) +
  scale_color_viridis(name = "Hair Coverage (%)",
                     option="plasma") +
  labs(x = "Longitude", y = "Latitude") +
  theme(panel.background = element_rect(fill = "white"))

ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group), fill = "gray", color = "gray", size = 0.1) +
  geom_point(data = df, aes(x = Longitude, y = Latitude, color = pilosity_percentage), size = 1) +
  scale_color_gradient(name = "Lightness",
                      low = "#180C3D",
                      high = "white") +
  labs(x = "Longitude", y = "Latitude") +
  theme(panel.background = element_rect(fill = "white"))

# blue continent
ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group), fill = "#B6C6AB", color = "#B6C6AB", size = 0.1) +
  geom_point(data = df, aes(x = Longitude, y = Latitude, color = pilosity_percentage), size = 1) +
  scale_color_viridis(name=NULL, option="plasma") +
  labs(x = "Longitude", y = "Latitude") +
  theme(panel.background = element_rect(fill = "white"))

ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group), fill = "#B6C6AB", color = "#B6C6AB", size = 0.1) +
  geom_point(data = df, aes(x = Longitude, y = Latitude, color = Brightness.Median), size = 1) +
  scale_color_gradient(name=NULL,
                       low = "#180C3D",
                       high = "white") +
  labs(x = "Longitude", y = "Latitude") +
  theme(panel.background = element_rect(fill = "white"))



#### Plot phylogenetic tree with trait values ####

# Calculate mean pilosity values for each Taxonomic_Info group
mean_pilosity <- aggregate(pilosity_percentage ~ Taxonomic_Info, data = matched_df, FUN = mean)

# Name the vector elements based on Taxonomic_Info
rownames(mean_pilosity) <- mean_pilosity$Taxonomic_Info

mean_pilosity$pilosity_percentage<- as.numeric(mean_pilosity$pilosity_percentage)
mean_pilosity_matrix <- as.matrix(mean_pilosity[, "pilosity_percentage", drop=FALSE])

# Calculate mean lightness values for each Taxonomic_Info group
mean_lightness <- aggregate(Brightness.Median ~ Taxonomic_Info, data = matched_df, FUN = mean)

# Name the vector elements based on Taxonomic_Info
rownames(mean_lightness) <- mean_lightness$Taxonomic_Info

mean_lightness$Brightness.Median<- as.numeric(mean_lightness$Brightness.Median)
mean_lightness_matrix <- as.matrix(mean_lightness[, "Brightness.Median", drop=FALSE])

# Plot tree
p <- ggtree(pruned_tree, branch.length='none', layout='circular')
#p <- ggtree(pruned_tree, branch.length='none', layout='circular') + geom_tiplab()

# Add the first heatmap
p1 <- gheatmap(p, mean_pilosity_matrix,
               offset=0.1,
               width=0.2,
               colnames=FALSE) + 
  scale_fill_viridis(name = "Hair Coverage (%)",
                    option="plasma")
# Add the first heatmap with the "spectral" color palette

# Add new scale for the second heatmap
p1 <- p1 + new_scale_fill()

# Add the second heatmap
p2 <- gheatmap(p1, mean_lightness_matrix,
               offset=5,  # Adjust the offset to avoid overlapping
               width=0.2,
               colnames=FALSE) + 
  scale_fill_continuous(name = "Lightness",
                        low = "#180C3D", high = "white",
                        na.value = "yellow")

# Display the combined plot
print(p2)


#### Plot trait values by family ####


ggplot(data=df, aes(x=Family, y=pilosity_percentage))+
  geom_boxplot(fill="#F26839")+
  labs(x = "Family", y = "Hair Coverage (%)") +  
  theme_classic() +
  theme(axis.text.x = element_text(color="black", size = 12, angle = 45, hjust = 1),axis.text.y = element_text(color="black", size = 12))


ggplot(data=df, aes(x=Family, y=Brightness.Median))+
  geom_boxplot(fill="#7D7EB3")+
  labs(x = "Family", y = "Lightness") +  
  theme_classic() +
  theme(axis.text.x = element_text(color="black", size = 12, angle = 45, hjust = 1),axis.text.y = element_text(color="black", size = 12))







