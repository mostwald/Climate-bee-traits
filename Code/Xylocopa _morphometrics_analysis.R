# Estimating body size in the large carpenter bees (Xylocopa)
# March 15, 2023

#### Load libraries and read in data frame ####
library(ggplot2)
library(lme4)
library(car)

xylo.morph.df <- read.csv("Xylocopa Morphometrics.csv")


#### Organize data ####


#add columns for log transformations
xylo.morph.df$log.dry.mass <- log(xylo.morph.df$Dry.Mass.without.pin..g.)
xylo.morph.df$log.dry.mass.mg <- log(xylo.morph.df$Dry.Mass.without.pin..mg.)
xylo.morph.df$log.itd <- log(xylo.morph.df$ITD..mm.)
xylo.morph.df$log.hw <- log(xylo.morph.df$Head.Width..mm.)
xylo.morph.df$log.cvl <- log(xylo.morph.df$Costal.Vein.Length..mm.)

#subset dataset to just include the columns we need
xylo.morph.df.subset <- subset(xylo.morph.df, select=c("Species","Sex","log.dry.mass","log.itd","log.hw","log.cvl"))
xylo.morph.df.subset.itd <- subset(xylo.morph.df, select=c("Species","Sex","log.dry.mass","log.itd"))
xylo.morph.df.subset.hw <- subset(xylo.morph.df, select=c("Species","Sex","log.dry.mass","log.hw"))
xylo.morph.df.subset.cvl <- subset(xylo.morph.df, select=c("Species","Sex","log.dry.mass","log.cvl"))

#subset data by sex
males <- subset(xylo.morph.df, Sex=="M")
females <- subset(xylo.morph.df, Sex=="F")





#### Model selection by stepwise regression ####
full.model.itd <- lm(formula = log.dry.mass ~ ., data = xylo.morph.df.subset.itd)
full.model.hw <- lm(formula = log.dry.mass ~ ., data = xylo.morph.df.subset.hw)
full.model.cvl <- lm(formula = log.dry.mass ~ ., data = xylo.morph.df.subset.cvl)

step.itd <- step(full.model.itd, scope = . ~ .^2, direction = 'both') # the ".^2" includes interactions among all terms.
summary(step.itd)
step.hw <- step(full.model.hw, scope = . ~ .^2, direction = 'both') #lowest AIC. Formula is mass ~ Species + Sex + HW + Species:Sex
summary(step.hw)
step.cvl <- step(full.model.cvl, scope = . ~ .^2, direction = 'both')
summary(step.cvl)


# Save selected models
selected.mod.itd <- lm(formula = log.dry.mass ~ Species + Sex + log.itd + Species:Sex + 
                         Sex:log.itd, data = xylo.morph.df.subset.itd)
summary(selected.mod.itd)
selected.mod.hw <- lm(formula = log.dry.mass ~ Species + Sex + log.hw + Species:Sex, 
                      data = xylo.morph.df.subset.hw)
selected.mod.cvl <- lm(formula = log.dry.mass ~ Species + Sex + log.cvl + Species:Sex + 
                         Sex:log.cvl, data = xylo.morph.df.subset.cvl)




#### Check assumptions of interspecific models (normality and homoscedasticity) ####

qqnorm(residuals(selected.mod.itd))
qqline(residuals(selected.mod.itd))
qqPlot(residuals(selected.mod.itd)) 
plot(fitted(selected.mod.itd),residuals(selected.mod.itd)) 

qqnorm(residuals(selected.mod.hw)) 
qqline(residuals(selected.mod.hw))
qqPlot(residuals(selected.mod.hw)) 
plot(fitted(selected.mod.hw),residuals(selected.mod.hw)) 

qqnorm(residuals(selected.mod.cvl)) 
qqline(residuals(selected.mod.cvl))
qqPlot(residuals(selected.mod.cvl)) 
plot(fitted(selected.mod.cvl),residuals(selected.mod.cvl)) 




#### How well do these size metrics predict intraspecific variation in body size? ####

# create subsets for each species that has > 8 individuals
x.tenuiscapa.females <- subset(females, Species =="Xylocopa tenuiscapa")
x.tranquebarica.females <- subset(females, Species =="Xylocopa tranquebarica")
x.tab.tab.females <- subset(females, Species =="Xylocopa tabaniformis tabaniformis")
x.tab.orp.females <- subset(females, Species =="Xylocopa tabaniformis orpifex")
x.frontalis.females <- subset(females, Species =="Xylocopa frontalis")
x.sonorina.females <- subset(females, Species =="Xylocopa sonorina")
x.inconstans.females <- subset(females, Species =="Xylocopa inconstans")
x.micans.females <- subset(females, Species =="Xylocopa micans")
x.californica.arizonensis.females <- subset(females, Species =="Xylocopa californica arizonensis")
x.virginica.females <- subset(females, Species =="Xylocopa virginica")
females.intraspecific.comparison <- rbind(x.virginica.females,x.tenuiscapa.females,x.tranquebarica.females,x.tab.tab.females,x.tab.orp.females,x.frontalis.females,x.sonorina.females,x.inconstans.females,x.micans.females,x.californica.arizonensis.females)
table(females$Species) # get sample size for each species

# OLS regression for each species for ITD
x.tenuiscapa.females.ols.itd <- lm(log.dry.mass ~ log.itd, data=x.tenuiscapa.females)
summary(x.tenuiscapa.females.ols.itd)
x.tranquebarica.females.ols.itd <- lm(log.dry.mass ~ log.itd, data=x.tranquebarica.females)
summary(x.tranquebarica.females.ols.itd)
x.tab.tab.females.ols.itd <- lm(log.dry.mass ~ log.itd, data=x.tab.tab.females)
summary(x.tab.tab.females.ols.itd)
x.tab.orp.females.ols.itd <- lm(log.dry.mass ~ log.itd, data=x.tab.orp.females)
summary(x.tab.orp.females.ols.itd)
x.frontalis.females.ols.itd <- lm(log.dry.mass ~ log.itd, data=x.frontalis.females)
summary(x.frontalis.females.ols.itd)
x.sonorina.females.ols.itd <- lm(log.dry.mass ~ log.itd, data=x.sonorina.females)
summary(x.sonorina.females.ols.itd)
x.inconstans.females.ols.itd <- lm(log.dry.mass ~ log.itd, data=x.inconstans.females)
summary(x.inconstans.females.ols.itd)
x.micans.females.ols.itd <- lm(log.dry.mass ~ log.itd, data=x.micans.females)
summary(x.micans.females.ols.itd)
x.californica.arizonensis.females.ols.itd <- lm(log.dry.mass ~ log.itd, data=x.californica.arizonensis.females)
summary(x.californica.arizonensis.females.ols.itd)
x.virginica.females.ols.itd <- lm(log.dry.mass ~ log.itd, data=x.virginica.females)
summary(x.virginica.females.ols.itd)

# OLS regression for each species for HW
x.tenuiscapa.females.ols.hw <- lm(log.dry.mass ~ log.hw, data=x.tenuiscapa.females)
summary(x.tenuiscapa.females.ols.hw)
x.tranquebarica.females.ols.hw <- lm(log.dry.mass ~ log.hw, data=x.tranquebarica.females)
summary(x.tranquebarica.females.ols.hw)
x.tab.tab.females.ols.hw <- lm(log.dry.mass ~ log.hw, data=x.tab.tab.females)
summary(x.tab.tab.females.ols.hw)
x.tab.orp.females.ols.hw <- lm(log.dry.mass ~ log.hw, data=x.tab.orp.females)
summary(x.tab.orp.females.ols.hw)
x.frontalis.females.ols.hw <- lm(log.dry.mass ~ log.hw, data=x.frontalis.females)
summary(x.frontalis.females.ols.hw)
x.sonorina.females.ols.hw <- lm(log.dry.mass ~ log.hw, data=x.sonorina.females)
summary(x.sonorina.females.ols.hw)
x.inconstans.females.ols.hw <- lm(log.dry.mass ~ log.hw, data=x.inconstans.females)
summary(x.inconstans.females.ols.hw)
x.micans.females.ols.hw <- lm(log.dry.mass ~ log.hw, data=x.micans.females)
summary(x.micans.females.ols.hw)
x.californica.arizonensis.females.ols <- lm(log.dry.mass ~ log.hw, data=x.californica.arizonensis.females)
summary(x.californica.arizonensis.females.ols)
x.virginica.females.ols.hw <- lm(log.dry.mass ~ log.hw, data=x.virginica.females)
summary(x.virginica.females.ols.hw)

# OLS regression for each species for CVL
x.tenuiscapa.females.ols.cvl <- lm(log.dry.mass ~ log.cvl, data=x.tenuiscapa.females)
summary(x.tenuiscapa.females.ols.cvl)
x.tranquebarica.females.ols.cvl <- lm(log.dry.mass ~ log.cvl, data=x.tranquebarica.females)
summary(x.tranquebarica.females.ols.cvl)
x.tab.tab.females.ols.cvl <- lm(log.dry.mass ~ log.cvl, data=x.tab.tab.females)
summary(x.tab.tab.females.ols.cvl)
x.tab.orp.females.ols.cvl <- lm(log.dry.mass ~ log.cvl, data=x.tab.orp.females)
summary(x.tab.orp.females.ols.cvl)
x.frontalis.females.ols.cvl <- lm(log.dry.mass ~ log.cvl, data=x.frontalis.females)
summary(x.frontalis.females.ols.cvl)
x.sonorina.females.ols.cvl <- lm(log.dry.mass ~ log.cvl, data=x.sonorina.females)
summary(x.sonorina.females.ols.cvl)
x.inconstans.females.ols.cvl <- lm(log.dry.mass ~ log.cvl, data=x.inconstans.females)
summary(x.inconstans.females.ols.cvl)
x.micans.females.ols.cvl <- lm(log.dry.mass ~ log.cvl, data=x.micans.females)
summary(x.micans.females.ols.cvl)
x.californica.arizonensis.females.ols.cvl <- lm(log.dry.mass ~ log.cvl, data=x.californica.arizonensis.females)
summary(x.californica.arizonensis.females.ols.cvl)
x.virginica.females.ols.cvl <- lm(log.dry.mass ~ log.cvl, data=x.virginica.females)
summary(x.virginica.females.ols.cvl)








#### Check assumptions of intraspecific models (normality and homoscedasticity) ####


# ITD
qqPlot(residuals(x.tenuiscapa.females.ols.itd)) 
qqPlot(residuals(x.tranquebarica.females.ols.itd)) 
qqPlot(residuals(x.tab.tab.females.ols.itd)) 
qqPlot(residuals(x.tab.orp.females.ols.itd)) 
qqPlot(residuals(x.frontalis.females.ols.itd))
qqPlot(residuals(x.sonorina.females.ols.itd)) # a couple outliers
qqPlot(residuals(x.micans.females.ols.itd)) 
qqPlot(residuals(x.inconstans.females.ols.itd)) 
qqPlot(residuals(x.californica.arizonensis.females.ols.itd)) 
qqPlot(residuals(x.virginica.females.ols.itd)) 

plot(fitted(x.tenuiscapa.females.ols.itd),residuals(x.tenuiscapa.females.ols.itd))
plot(fitted(x.tranquebarica.females.ols.itd),residuals(x.tranquebarica.females.ols.itd))
plot(fitted(x.tab.tab.females.ols.itd),residuals(x.tab.tab.females.ols.itd))
plot(fitted(x.tab.orp.females.ols.itd),residuals(x.tab.orp.females.ols.itd))
plot(fitted(x.frontalis.females.ols.itd),residuals(x.frontalis.females.ols.itd))
plot(fitted(x.sonorina.females.ols.itd),residuals(x.sonorina.females.ols.itd))
plot(fitted(x.micans.females.ols.itd),residuals(x.micans.females.ols.itd))
plot(fitted(x.inconstans.females.ols.itd),residuals(x.inconstans.females.ols.itd))
plot(fitted(x.californica.arizonensis.females.ols.itd),residuals(x.californica.arizonensis.females.ols.itd))
plot(fitted(x.virginica.females.ols.itd),residuals(x.virginica.females.ols.itd))


# HW
qqPlot(residuals(x.tenuiscapa.females.ols.hw)) 
qqPlot(residuals(x.tranquebarica.females.ols.hw)) 
qqPlot(residuals(x.tab.tab.females.ols.hw)) 
qqPlot(residuals(x.tab.orp.females.ols.hw)) # a couple outliers
qqPlot(residuals(x.frontalis.females.ols.hw))
qqPlot(residuals(x.sonorina.females.ols.hw)) # a couple outliers
qqPlot(residuals(x.micans.females.ols.hw)) 
qqPlot(residuals(x.inconstans.females.ols.hw)) 
qqPlot(residuals(x.californica.arizonensis.females.ols.hw)) 
qqPlot(residuals(x.virginica.females.ols.hw)) 

plot(fitted(x.tenuiscapa.females.ols.hw),residuals(x.tenuiscapa.females.ols.hw))
plot(fitted(x.tranquebarica.females.ols.hw),residuals(x.tranquebarica.females.ols.hw))
plot(fitted(x.tab.tab.females.ols.hw),residuals(x.tab.tab.females.ols.hw))
plot(fitted(x.tab.orp.females.ols.hw),residuals(x.tab.orp.females.ols.hw))
plot(fitted(x.frontalis.females.ols.hw),residuals(x.frontalis.females.ols.hw))
plot(fitted(x.sonorina.females.ols.hw),residuals(x.sonorina.females.ols.hw))
plot(fitted(x.micans.females.ols.hw),residuals(x.micans.females.ols.hw))
plot(fitted(x.inconstans.females.ols.hw),residuals(x.inconstans.females.ols.hw))
plot(fitted(x.californica.arizonensis.females.ols.hw),residuals(x.californica.arizonensis.females.ols.hw))
plot(fitted(x.virginica.females.ols.hw),residuals(x.virginica.females.ols.hw))

# CVL
qqPlot(residuals(x.tenuiscapa.females.ols.cvl)) 
qqPlot(residuals(x.tranquebarica.females.ols.cvl)) 
qqPlot(residuals(x.tab.tab.females.ols.cvl)) 
qqPlot(residuals(x.tab.orp.females.ols.cvl)) 
qqPlot(residuals(x.frontalis.females.ols.cvl))
qqPlot(residuals(x.sonorina.females.ols.cvl)) # outliers
qqPlot(residuals(x.micans.females.ols.cvl)) 
qqPlot(residuals(x.inconstans.females.ols.cvl)) 
qqPlot(residuals(x.californica.arizonensis.females.ols.cvl)) 
qqPlot(residuals(x.virginica.females.ols.cvl)) 

plot(fitted(x.tenuiscapa.females.ols.cvl),residuals(x.tenuiscapa.females.ols.cvl))
plot(fitted(x.tranquebarica.females.ols.cvl),residuals(x.tranquebarica.females.ols.cvl))
plot(fitted(x.tab.tab.females.ols.cvl),residuals(x.tab.tab.females.ols.cvl))
plot(fitted(x.tab.orp.females.ols.cvl),residuals(x.tab.orp.females.ols.cvl))
plot(fitted(x.frontalis.females.ols.cvl),residuals(x.frontalis.females.ols.cvl))
plot(fitted(x.sonorina.females.ols.cvl),residuals(x.sonorina.females.ols.cvl))
plot(fitted(x.micans.females.ols.cvl),residuals(x.micans.females.ols.cvl))
plot(fitted(x.inconstans.females.ols.cvl),residuals(x.inconstans.females.ols.cvl))
plot(fitted(x.californica.arizonensis.females.ols.cvl),residuals(x.californica.arizonensis.females.ols.cvl))
plot(fitted(x.virginica.females.ols.cvl),residuals(x.virginica.females.ols.cvl))










#### Plot data ####

# Plot log dry mass vs. log ITD for all bees
ggplot(xylo.morph.df, aes(x=log.dry.mass.mg, y=log.itd, color=Species, shape=Sex)) +
  geom_point() +
  theme_classic() +
  xlab("ln Mass (mg)") +
  ylab("ln ITD (mm)") +
  theme(axis.text=element_text(size=12, color="black"),axis.title=element_text(size=12, color="black"), legend.position = "none")

# Plot log dry mass vs. log HW for all bees
ggplot(xylo.morph.df, aes(x=(log.dry.mass.mg), y=(log.hw), color=Species, shape=Sex)) +
  geom_point() +
  theme_classic() +
  xlab("ln Mass (mg)") +
  ylab("ln Head Width (mm)") +
  theme(axis.text=element_text(size=12, color="black"),axis.title=element_text(size=12, color="black"),legend.position="none") 

# Plot log dry mass vs. log CVL for all bees
ggplot(xylo.morph.df, aes(x=(log.dry.mass.mg), y=(log.cvl), color=Species, shape=Sex)) +
  geom_point() +
  theme_classic() +
  xlab("ln Mass (mg)") +
  ylab("ln Costal Vein Length (mm)") +
  theme(axis.text=element_text(size=12, color="black"),axis.title=element_text(size=12, color="black"),legend.position="none")

# find out what colors were used above:
library(scales)
for(i in 1:14){
  print(hue_pal()(i))
}
# x. cal ariz = "#F8766D"
# x. cal cal = "#E38900"
# x. cal diamesa = "#C49A00"
# x. frontalis = "#99A800"
# x. inconstans = "#53B400"
# x. latipes = "#00BC56"
# x. micans = "#00C094"
# x. ruficornis = "#00BFC4"
# X. sonorina = "#00B6EB"
# x. tab. orp. = "#06A4FF"
# x. tab. tab = "#A58AFF"
# x. tenuiscapa = "#DF70F8"
# x. tranquebarica = "#FB61D7"
# x. virginica ="#FF66A8"



# Plot log mass vs log ITD for females only, with regression lines
ggplot(females.intraspecific.comparison, aes(x=log.dry.mass, y=log.itd, color=Species)) +
  geom_point() +
  theme_classic() +
  xlab("ln Mass (mg)") +
  ylab("ln ITD (mm)") +
  geom_smooth(method="lm") +  # gray areas are 95% CI
  theme(axis.text=element_text(size=12, color="black"),axis.title=element_text(size=12, color="black"), legend.position="none") +
  scale_color_manual(values=c("#F8766D","#99A800","#53B400", "#00C094", "#00B6EB", "#06A4FF", "#A58AFF", "#DF70F8", "#FB61D7", "#FF66A8"))  +
  ylim(1.45,2.35)

# Plot log mass vs log HW for females only, with regression lines
ggplot(females.intraspecific.comparison, aes(x=log.dry.mass, y=log.hw, color=Species)) +
  geom_point() +
  theme_classic() +
  xlab("ln Mass (mg)") +
  ylab("ln Head Width (mm)") +
  geom_smooth(method="lm") +  # gray areas are 95% CI
  theme(axis.text=element_text(size=12, color="black"),axis.title=element_text(size=12, color="black"), legend.position = "none") +
  scale_color_manual(values=c("#F8766D","#99A800","#53B400", "#00C094", "#00B6EB", "#06A4FF", "#A58AFF", "#DF70F8", "#FB61D7", "#FF66A8")) +
  ylim(1.45,2.35)

# Plot log mass vs log CVL for females only, with regression lines
ggplot(females.intraspecific.comparison, aes(x=log.dry.mass, y=log.cvl, color=Species)) +
  geom_point() +
  theme_classic() +
  xlab("ln Mass (mg)") +
  ylab("ln Costal Vein Length (mm)") +
  geom_smooth(method="lm") +  # gray areas are 95% CI
  theme(axis.text=element_text(size=12, color="black"),axis.title=element_text(size=12, color="black"), legend.position = "none") +
  scale_color_manual(values=c("#F8766D","#99A800","#53B400", "#00C094", "#00B6EB", "#06A4FF", "#A58AFF", "#DF70F8", "#FB61D7", "#FF66A8"))  +
  ylim(1.45,2.35)





#### Analyze 3D model data ####

species <- c("X. virginica", 
             "X. micans",
             "X. californica arizonensis",
             "X. sonorina",
             "X. micans",
             "X. tenuiscapa",
             "X. tabaniformis orpifex",
             "X. ruficornis",
             "X. inconstans",
             "X. frontalis",
             "X. tabniformis tabaniformis",
             "X. latipes",
             "X. tranquebarica")
catalog_number <- c("LACM ENT 602591",
                    "LACM ENT_602478",
                    "LACM ENT_602461",
                    "UCSB-IZC00012185",
                    "LACM ENT 602464",
                    "LACM ENT 602432",
                    "UCSB-IZC00014711",
                    "LACM ENT 602427",
                    "LACM ENT 602500",
                    "LACM ENT 602402",
                    "LACM ENT 602423",
                    "LACM ENT 601466",
                    "LACM ENT 601584")
log_dry_mass <- c(5.185708,
                  5.015291,
                  5.757955,
                  5.84846,
                  5.098035,
                  6.47497,
                  4.94663,
                  6.025141,
                  6.022721,
                  6.515749,
                  4.872905,
                  6.814214,
                  6.82622)
volume <- c(875.7,
            477.947,
            699.462,
            1575,
            659.708,
            2419,
            530.583,
            2385.6,
            663.826,
            3118,
            674.1,
            4614.4,
            2262.8
)
surface_area <- c(1126.3,
                  678.403,
                  811.785,
                  1616.8,
                  780.728,
                  1855,
                  759.653,
                  1962.7,
                  788.019,
                  2148.9,
                  922.3,
                  2275,
                  1843.5)
threeD_matrix <- cbind(species, catalog_number, log_dry_mass, volume, surface_area)
threeD_df <- as.data.frame(threeD_matrix)
threeD_df$volume <- as.numeric(threeD_df$volume)
threeD_df$log_dry_mass <- as.numeric(threeD_df$log_dry_mass)
threeD_df$surface_area <- as.numeric(threeD_df$surface_area)
threeD_df$log_volume <- log(threeD_df$volume)
threeD_df$log_surface_area <- log(threeD_df$surface_area)

# Plot log mass vs log volume 
ggplot(threeD_df, aes(x=log_dry_mass, y=log_volume, color=species)) +
  geom_point() +
  theme_classic() +
  xlab("ln Mass (mg)") +
  ylab(expression(paste("ln volume (mm"^3*")"))) +
  geom_smooth(data = threeD_df, aes(x = log_dry_mass, y = log_volume), method="lm", linetype = "solid", color = "black") +
  theme(axis.text=element_text(size=12, color="black"),axis.title=element_text(size=12, color="black")) +
  scale_color_manual(values=c("#F8766D","#99A800","#53B400", "#00BC56", "#00C094", "#00BFC4", "#00B6EB","#06A4FF","#A58AFF","#DF70F8","#FB61D7","#FF66A8"))
ggplot(threeD_df, aes(x=log_dry_mass, y=log_surface_area, color=species)) +
  geom_point() +
  theme_classic() +
  xlab("ln Mass (mg)") +
  ylab(expression(paste("ln surface area (mm"^2*")"))) +
  geom_smooth(data = threeD_df, aes(x = log_dry_mass, y = log_volume), method="lm", linetype = "solid", color = "black") +
  theme(axis.text=element_text(size=12, color="black"),axis.title=element_text(size=12, color="black"),legend.text=element_text(face="italic")) +
  scale_color_manual(values=c("#F8766D","#99A800","#53B400", "#00BC56", "#00C094", "#00BFC4", "#00B6EB","#06A4FF","#A58AFF","#DF70F8","#FB61D7","#FF66A8"))

# OLS regressions for surface area and volume
surface.area.ols <- lm(log_dry_mass ~ log_surface_area, data=threeD_df)
volume.ols <- lm(log_dry_mass ~ log_volume, data=threeD_df)
summary(surface.area.ols)
summary(volume.ols)


qqPlot(residuals(surface.area.ols)) 
plot(fitted(surface.area.ols),residuals(surface.area.ols))
qqPlot(residuals(volume.ols)) 
plot(fitted(volume.ols),residuals(volume.ols))
