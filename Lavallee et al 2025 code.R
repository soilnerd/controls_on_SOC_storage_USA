### Code for Lavallee et al. 2025 ###
# 0) Read data ----
library(tidyverse)
setwd("/Users/jocelynlavallee/Dropbox/3\ Other\ Projects/NRCS\ Soil\ Inventory/Manuscript")
load("Lavallee_et_al_2025_data.RData")

# 0) Characteristics of the dataset ----
# Land cover - note all samples with an Ap horizon were manually classified as "cropland"
summary_meas <- dat %>%
  group_by(fraction_quant_method) %>%
  summarize(n())

summary_cover <- dat %>%
  group_by(land_cover_simple, horizon_ap) %>%
  summarize(n())

summary_NPP_range <- dat %>%
  group_by(AI_factor) %>%
  summarize(min(NPP), max(NPP)) %>%
  mutate(range = `max(NPP)` - `min(NPP)`)


# 1) Regression of SOC and NPP vs. Aridity Index ----
library(segmented)
dat_npp_soc_lm <- dat %>%
  dplyr::select(WS_mgpg_OC, NPP, AI) %>%
  drop_na() 

#linear model
#R2 = 0.321, Adjusted R2 = 0.321
lm_npp_ai <- lm(NPP ~ AI, data = dat_npp_soc_lm)
summary(lm_npp_ai)

# Suggesting two breakpoints gives 0.42 and 1.37
# R2 = 0.536, Adjusted R2 = 0.534
seg_npp_ai_free <- segmented(lm_npp_ai, 
                             seg.Z = ~ AI, 
                             npsi = 2,
                             psi = c(0.5, 0.8))
# display the summary
summary(seg_npp_ai_free)

# Providing fixed breakpoint at AI = 0.65
# R2 = 0.5, Adjusted R2 = 0.50
seg_npp_ai_65 <- segmented(lm_npp_ai, 
                           seg.Z = ~ AI, 
                           psi = 0.65,
                           control = seg.control(it.max = 0)
                           )
summary(seg_npp_ai_65)
# slope = 0.88373 + -.79606 = 0.08767

# Compare to a penalized cubic regression spline GAM (default when plotted with smooth function in ggplot)
# in mgcv::gam() is used with formula = y ~ s(x, bs = "cs") with method = "REML"
# R2 = 0.534
library(mgcv)
gam_npp_ai <- mgcv::gam(NPP ~ s(AI, bs = "cs"), data = dat_npp_soc_lm, method = "REML")
summary(gam_npp_ai)


# Get the fitted data
fitted_npp_free <- fitted(seg_npp_ai_free)
segmodel_npp_free <- data.frame(AI = dat_npp_soc_lm$AI, NPP = fitted_npp_free)
fitted_npp_65 <- fitted(seg_npp_ai_65)
segmodel_npp_65 <- data.frame(AI = dat_npp_soc_lm$AI, NPP = fitted_npp_65)

# plot the fitted NPP models
ggplot(segmodel_npp_65, aes(x = AI, y = NPP)) + 
  geom_line() +
  geom_line(data = segmodel_npp_free,  aes(x = AI, y = NPP))


# Analysis for SOC
# Linear model
# Multiple R2 = 0.120, Adjusted R2 = 0.119
lm_soc_ai <- lm(WS_mgpg_OC~ AI, data = dat_npp_soc_lm)
summary(lm_soc_ai)

# Suggesting two breakpoints gives 0.45 and 0.79
# R2 = 0.158, Adjusted R2 = 154.
seg_soc_ai_free <- segmented(lm_soc_ai, 
                             seg.Z = ~ AI, 
                             npsi = 2,
                             psi = c(0.5, 0.8))
summary(seg_soc_ai_free)

# Forcing the breakpoint at AI = 0.65
# R2 = 0.125, Adjusted R2 = 0.123
seg_soc_ai_65 <- segmented(lm_soc_ai, 
                           seg.Z = ~ AI, 
                           psi = 0.65,
                           fixed.psi = 0.65,
                           control = seg.control(it.max = 0))
summary(seg_soc_ai_65)
# slope = 19.088 + -9.246 = 9.842


# Compare to a penalized cubic regression spline GAM (default when plotted with smooth function in ggplot)
# in mgcv::gam() is used with formula = y ~ s(x, bs = "cs") with method = "REML"
# R2 = 0.169
gam_soc_ai <- mgcv::gam(WS_mgpg_OC ~ s(AI, bs = "cs"), data = dat_npp_soc_lm, method = "REML")
summary(gam_soc_ai)

# Get the fitted data
fitted_soc_free <- fitted(seg_soc_ai_free)
segmodel_soc_free <- data.frame(AI = dat_npp_soc_lm$AI, SOC = fitted_soc_free)
fitted_soc_65 <- fitted(seg_soc_ai_65)
segmodel_soc_65 <- data.frame(AI = dat_npp_soc_lm$AI, SOC = fitted_soc_65)

# plot the fitted SOC models
ggplot(segmodel_soc_65, aes(x = AI, y = SOC)) + 
  geom_line() +
  geom_line(data = segmodel_soc_free,  aes(x = AI, y = SOC))


## Fig. 1 ----
library(ggpmisc)
Fig_1A <-  
  ggplot(data = dat,
         aes(x = AI, y = NPP)) + 
  geom_point(alpha = .4, col = "#49beaa") +
  geom_vline(xintercept = 0.65, color = "#9381ff", linetype = 2) +
  geom_vline(xintercept = 0.42, color = "#eeb868", linetype = 2) +
  geom_vline(xintercept = 1.37, color = "#eeb868", linetype = 2) +
  # plot the fitted models
  geom_smooth(aes(y = NPP), color = "#595959", method = "gam", alpha = 0.2, linewidth = 0.5) +
  geom_line(data = segmodel_65, aes(x = AI, y = NPP), col = "#9381ff", linewidth = 0.7) + ##185e40
  geom_line(data = segmodel_free,  aes(x = AI, y = NPP), col = "#eeb868", linewidth = 0.8) + ##0c2f20
  scale_y_continuous(expand = c(0, 0), name = expression(paste("NPP (kg C m"^-2, ")")),
                     limits = c(0, 2.3)) +
  scale_x_continuous(breaks = c(0, 0.42, 0.65, 1.37, 1, 2, 3, 4),
                     labels = c("0", "0.42","0.65", "1.37", "1", "2", "3", "4"),
                     expand = c(0, 0), limits = c(0,4), name = "Aridity Index") + # Omits 3 points > 4.
  annotate("text", x = 3.5, y = 1.5, label = "R^2 == 0.54", col = "#eeb868", parse=TRUE) +
  annotate("text", x = 3.5, y = 1.3, label = "R^2 == 0.50", col = "#9381ff", parse=TRUE) +
  theme_classic()

Fig_1B <- 
  ggplot(data = dat,
         aes(x = AI)) +
  geom_point(aes(y = WS_mgpg_OC), color = "#eeb868", 
             alpha = .5) +
  geom_vline(xintercept = 0.65, color = "#ba2c73", linetype = 2) +
  geom_vline(xintercept = 0.45, color = "#49beaa", linetype = 2) +
  geom_vline(xintercept = 0.79, color = "#49beaa", linetype = 2) +
  #plot the fitted models
  geom_smooth(aes(y = WS_mgpg_OC), color = "#595959", method = "gam", alpha = 0.2, linewidth = 0.5) +
  geom_line(data = segmodel_soc_65, aes(x = AI, y = SOC), col = "#ba2c73", linewidth = 0.7) + ##185e40
  geom_line(data = segmodel_soc_free,  aes(x = AI, y = SOC), col = "#49beaa", linewidth = 0.8) + ##0c2f20

  scale_y_continuous(expand = c(0, 0), name = expression(paste("SOC (g C kg"^-1, "soil)"))) +
  scale_x_continuous(breaks = c(0, 0.45, 0.65, 0.79, 1, 2, 3, 4),
                     labels = c("0", "0.45", "0.65", "0.79", "1", "2", "3", "4"),
                     expand = c(0, 0), limits = c(0,4), name = "Aridity Index") + # Omits 4 points > 4
  annotate("text", x = 3.5, y = 30, label = "R^2 == 0.16", col = "#49beaa", parse=TRUE) +
  annotate("text", x = 3.5, y = 23, label = "R^2 == 0.13", col = "#ba2c73", parse=TRUE) +
  theme_classic()

Fig_1C <- 
  ggplot(aes(x = AI_factor), dat = dat) +
  geom_point(aes(y = log(SOC_NPP_Ratio)), color = "#6a9abb", alpha = 0.4, position = "jitter",
             size = 1) +
  geom_boxplot(aes(y = log(SOC_NPP_Ratio)), color = "#064672", alpha = 0) +
  labs(x = "Aridity Index", y =  
         expression(paste("ln(SOC:NPP ratio)"))) +
  theme_classic() +
  theme(axis.title.y = element_text(color = "#064672")) +
  scale_y_continuous(expand = c(0, 0), limits = c(1, 6.5)) + # Omits 9 points outside range for readability
  annotate("text", x = c(1, 2), y=c(rep(6.3, 2)), label = c("a", "b"))

cat("Fig 1C  N = ", nrow(filter(dat, log(SOC_NPP_Ratio) > 1))) # no points > 6.5

library(patchwork)
layout <- "
AAA
BBB
CCC
"
p1 <- Fig_1A + Fig_1B + Fig_1C +
  plot_layout(design = layout) +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 16))
pdf('Fig_1.pdf', 
    width = 4, 
    height= 10, 
    useDingbats=FALSE)
p1
dev.off()


# 2) Regressions of NPP, POC, and MAOC vs. clay+silt, al_fe, ca_mg, MCI ----
## Fig. 2 ----
# NPP plots
plot_NPP_cs <- ggplot(dat, aes(y = NPP, x = claysilt)) +
  geom_point(aes(color = ph_h2o), size = 0.3) +
  scale_color_gradient(name = "pH", low = "orange", high = "blue") +
  stat_poly_line(linewidth = 0.5, color = "black") +
  stat_correlation(mapping = use_label("R"),
                   label.x = "right",
                   label.y = .95)+
  stat_correlation(mapping = use_label("R2", "P"),
                   label.x = "right",
                   label.y = .9) +
  theme_classic() +
  facet_grid(cols = vars(AI_factor), scales = "free") +
  ggtitle("Silt + Clay") +
  theme(plot.title = element_text(hjust = 0.5, size = 12)) +  
  ylab(expression(paste("NPP (kg C m"^-2, ")"))) +
  xlab(element_blank())

plot_NPP_ca_mg <- ggplot(dat, aes(y = NPP, x = ca_mg)) +
  geom_point(aes(color = ph_h2o), size = 0.3) +
  scale_color_gradient(name = "pH", low = "orange", high = "blue") +
  stat_poly_line(linewidth = 0.5, color = "black") +
  stat_correlation(mapping = use_label("R"),
                   label.x = "right",
                   label.y = .95)+
  stat_correlation(mapping = use_label("R2", "P"),
                   label.x = "right",
                   label.y = .9) +
  theme_classic() +
  facet_grid(cols = vars(AI_factor), scales = "free") +
  ggtitle(expression(Ca[ex] + Mg[ex])) +
  theme(plot.title = element_text(hjust = 0.5, size = 12)) +
  ylab(element_blank()) +
  xlab(element_blank())

plot_NPP_alFe <- ggplot(dat, aes(y = NPP, x = al_fe)) +
  geom_point(aes(color = ph_h2o), size = 0.3) +
  scale_color_gradient(name = "pH", low = "orange", high = "blue") +
  stat_poly_line(linewidth = 0.5, color = "black") +
  stat_correlation(mapping = use_label("R"),
                   label.x = "right",
                   label.y = .95)+
  stat_correlation(mapping = use_label("R2", "P"),
                   label.x = "right",
                   label.y = .9) +
  theme_classic() +
  facet_grid(cols = vars(AI_factor), scales = "free") +
  ggtitle(expression(paste(Al[o] + "[1/2]",Fe[o]))) +
  theme(plot.title = element_text(hjust = 0.5, size = 12)) +
  ylab(element_blank()) +
  xlab(element_blank())

plot_NPP_mci <- ggplot(dat, aes(y = NPP, x = MCI)) +
  geom_point(aes(color = ph_h2o), size = 0.3) +
  scale_color_gradient(name = "pH", low = "orange", high = "blue") +
  stat_poly_line(linewidth = 0.5, color = "black") +
  stat_correlation(mapping = use_label("R"),
                   label.x = "right",
                   label.y = .95)+
  stat_correlation(mapping = use_label("R2", "P"),
                   label.x = "right",
                   label.y = .9) +
  theme_classic() +
  facet_grid(cols = vars(AI_factor), scales = "free") +
  ggtitle("MCI") +
  theme(plot.title = element_text(hjust = 0.5, size = 12)) +
  ylab(element_blank()) +
  xlab(element_blank())

# POC plots
plot_POC_cs <- ggplot(dat, aes(y = POM_mgpg_OC, x = claysilt)) +
  geom_point(aes(color = ph_h2o), size = 0.3) +
  scale_color_gradient(name = "pH", low = "orange", high = "blue") +
  stat_poly_line(linewidth = 0.5, color = "black") +
  stat_correlation(mapping = use_label("R"),
                   label.x = "right",
                   label.y = .95)+
  stat_correlation(mapping = use_label("R2", "P"),
                   label.x = "right",
                   label.y = .9) +
  theme_classic() +
  facet_grid(cols = vars(AI_factor), scales = "free") +
  ylab(expression(paste("POC (g C kg"^-1, "soil)"))) +
  xlab(element_blank())

plot_POC_ca_mg <- ggplot(dat, aes(y = POM_mgpg_OC, x = ca_mg)) +
  geom_point(aes(color = ph_h2o), size = 0.3) +
  scale_color_gradient(name = "pH", low = "orange", high = "blue") +
  stat_poly_line(linewidth = 0.5, color = "black") +
  stat_correlation(mapping = use_label("R"),
                   label.x = "right",
                   label.y = .95)+
  stat_correlation(mapping = use_label("R2", "P"),
                   label.x = "right",
                   label.y = .9) +
  theme_classic() +
  facet_grid(cols = vars(AI_factor), scales = "free") +
  ylab(element_blank()) +
  xlab(element_blank())

plot_POC_alFe <- ggplot(dat, aes(y = POM_mgpg_OC, x = al_fe)) +
  geom_point(aes(color = ph_h2o), size = 0.3) +
  scale_color_gradient(name = "pH", low = "orange", high = "blue") +
  stat_poly_line(linewidth = 0.5, color = "black") +
  stat_correlation(mapping = use_label("R"),
                   label.x = "right",
                   label.y = .95)+
  stat_correlation(mapping = use_label("R2", "P"),
                   label.x = "right",
                   label.y = .9) +
  theme_classic() +
  facet_grid(cols = vars(AI_factor), scales = "free") +
  ylab(element_blank()) +
  xlab(element_blank())

plot_POC_mci <- ggplot(dat, aes(y = POM_mgpg_OC, x = MCI)) +
  geom_point(aes(color = ph_h2o), size = 0.3) +
  scale_color_gradient(name = "pH", low = "orange", high = "blue") +
  stat_poly_line(linewidth = 0.5, color = "black") +
  stat_correlation(mapping = use_label("R"),
                   label.x = "right",
                   label.y = .95)+
  stat_correlation(mapping = use_label("R2", "P"),
                   label.x = "right",
                   label.y = .9) +
  theme_classic() +
  facet_grid(cols = vars(AI_factor), scales = "free") +
  theme(plot.title = element_text(hjust = 0.5, size = 12)) +
  ylab(element_blank()) +
  xlab(element_blank())


#MAOC plots (bottom)
plot_MAOC_cs <- ggplot(dat, aes(y = MAOM_mgpg_OC, x = claysilt)) +
  geom_point(aes(color = ph_h2o), size = 0.3) +
  scale_color_gradient(name = "pH", low = "orange", high = "blue") +
  stat_poly_line(linewidth = 0.5, color = "black") +
  stat_correlation(mapping = use_label("R"),
                   label.x = "right",
                   label.y = .95)+  
  stat_correlation(mapping = use_label("R2", "P"),
                   label.x = "right",
                   label.y = .9) +
  theme_classic() +
  facet_grid(cols = vars(AI_factor), scales = "free") +
  ylab(expression(paste("MAOC (g C kg"^-1, "soil)"))) +
  xlab("(% dry soil mass)")

plot_MAOC_ca_mg <- ggplot(dat, aes(y = MAOM_mgpg_OC, x = ca_mg)) +
  geom_point(aes(color = ph_h2o), size = 0.3) +
  scale_color_gradient(name = "pH", low = "orange", high = "blue") +
  stat_poly_line(linewidth = 0.5, color = "black") +
  stat_correlation(mapping = use_label("R"),
                   label.x = "right",
                   label.y = .95)+
  stat_correlation(mapping = use_label("R2", "P"),
                   label.x = "right",
                   label.y = .9) +
  theme_classic() +
  facet_grid(cols = vars(AI_factor), scales = "free") +
  ylab(element_blank()) +
  xlab(expression(paste("cmol kg"^-1, "soil")))

plot_MAOC_alFe <- ggplot(dat, aes(y = MAOM_mgpg_OC, x = al_fe)) +
  geom_point(aes(color = ph_h2o), size = 0.3) +
  scale_color_gradient(name = "pH", low = "orange", high = "blue") +
  stat_poly_line(linewidth = 0.5, color = "black") +
  stat_correlation(mapping = use_label("R"),
                   label.x = "right",
                   label.y = .95)+
  stat_correlation(mapping = use_label("R2", "P"),
                   label.x = "right",
                   label.y = .9) +
  theme_classic() +
  facet_grid(cols = vars(AI_factor), scales = "free") +
  ylab(element_blank()) +
  xlab(expression(paste("g kg"^-1, "soil")))

plot_MAOC_mci <- ggplot(dat, aes(y = MAOM_mgpg_OC, x = MCI)) +
  geom_point(aes(color = ph_h2o), size = 0.3) +
  scale_color_gradient(name = "pH", low = "orange", high = "blue") +
  stat_poly_line(linewidth = 0.5, color = "black") +
  stat_correlation(mapping = use_label("R"),
                   label.x = "right",
                   label.y = .95)+
  stat_correlation(mapping = use_label("R2", "P"),
                   label.x = "right",
                   label.y = .9) +
  theme_classic() +
  facet_grid(cols = vars(AI_factor), scales = "free") +
  theme(plot.title = element_text(hjust = 0.5, size = 12)) +
  ylab(element_blank()) +
  xlab("Unitless")

#library(patchwork)
p2 <-   (plot_NPP_cs | plot_NPP_ca_mg | plot_NPP_alFe | plot_NPP_mci) /
  (plot_POC_cs | plot_POC_ca_mg | plot_POC_alFe | plot_POC_mci) /
  (plot_MAOC_cs | plot_MAOC_ca_mg | plot_MAOC_alFe | plot_MAOC_mci) +
  plot_layout(guides = 'collect')

pdf('Fig_2.pdf', 
    width = 16, 
    height= 8, 
    useDingbats=FALSE)
p2
dev.off()

# 3) Structural Equation Model ----
library(lavaan)
model_SEM <- '
POM_mgpg_OC ~ 
  c("sc.poc1", "sc.poc2") * claysilt + #dont constrain
  c("ph.poc1", "ph.poc1") * ph_h2o + # NS p = 0.587
  c("af.poc1", "af.poc2") * al_fe + # NS p = 0.05506, AIC lower by 2
  c("ca_mg.poc1", "ca_mg.poc2") * ca_mg + # p = 0.0017
  c("et.poc1", "et.poc1") * ET + # NS p = 0.2651
  c("map.poc1", "map.poc2") * prcp_mean +  # dont constrain, NS in arid
  c("npp.poc1", "npp.poc2") * NPP + # dont constrain, NS in humid
  c("cult.poc1", "cult.poc1") * horizon_ap + # NS p = 0.3324
  c("depth.poc1", "depth.poc1") * hzn_bot # NS p = 0.78
MAOM_mgpg_OC ~ 
  c("sc.maoc1", "sc.maoc1") * claysilt + # NS p = 0.95
  c("ph.maoc1", "ph.maoc2") * ph_h2o + # p = 0.011
  c("af.maoc1", "af.maoc1") * al_fe + # p = 0.402
  c("ca_mg.maoc1", "ca_mg.maoc2") * ca_mg + # p = 1.269e-05
  c("et.maoc1", "et.maoc2") * ET + # p = 0.04
  c("map.maoc1", "map.maoc2") * prcp_mean +  #dont constrain, NS in humid
  c("npp.maoc1", "npp.maoc2") * NPP + #dont constrain, NS in humid
  c("cult.maoc1", "cult.maoc1") * horizon_ap + # NS p = 0.6
  c("depth.maoc1", "depth.maoc1") * hzn_bot # NS p = 0.63
ph_h2o ~ 
  c("et.ph1", "et.ph2") * ET + # p = 0.001
  c("map.ph1", "map.ph2") * prcp_mean + # p = 8.98e-09
  c("sc.ph1", "sc.ph2") * claysilt + #dont constrain, NS in humid
  c("af.ph1", "af.ph2") * al_fe + #dont constrain, NS in arid
  c("ca_mg.ph1", "ca_mg.ph2") * ca_mg + # p = 0.01
  c("cult.ph1", "cult.ph2") * horizon_ap + #dont constrain, NS in arid
  c("depth.ph1", "depth.ph2") * hzn_bot #dont constrain, NS in arid
NPP ~ 
  c("ph.npp1", "ph.npp1") * ph_h2o + # NS p = 0.8
  c("et.npp1", "et.npp2") * ET + # p = 3.103e-07
  c("map.npp1", "map.npp2") * prcp_mean + # p = 3.055e-10
  c("af.npp1", "af.npp2") * al_fe + # dont constrain, NS arid
  c("ca_mg.npp1", "ca_mg.npp2") * ca_mg + # dont constrain, NS humid
  c("cult.npp1", "cult.npp2") * horizon_ap + # dont constrain, NS humid
  c("depth.npp1", "depth.npp2") * hzn_bot #dont constrain, NS arid
             
    # Differences, paths of interest
    # see https://groups.google.com/g/lavaan/c/I2Svn1gfTTU/m/Ma5Oouz2AAAJ
    cult.pom.v.maom := cult.poc1 - cult.maoc1 
    hzn.pom.v.maom := depth.poc1 - depth.maoc1
    ph.pom.v.maom1 := ph.poc1 - ph.maoc1
    ph.pom.v.maom2 := ph.poc1 - ph.maoc2
    af.pom.v.maom1 := af.poc1 - af.maoc1
    af.pom.v.maom2 := af.poc2 - af.maoc1
    cm.pom.v.maom1 := ca_mg.poc1 - ca_mg.maoc1
    cm.pom.v.maom2 := ca_mg.poc2 - ca_mg.maoc2    
    et.pom.v.maom1 := et.poc1 - et.maoc1
    et.pom.v.maom2 := et.poc1 - et.maoc2
    map.pom.v.maom1 := map.poc1 - map.maoc1
    map.pom.v.maom2 := map.poc2 - map.maoc2
    npp.pom.v.maom1 := npp.poc1 - npp.maoc1
    npp.pom.v.maom2 := npp.poc2 - npp.maoc2
    clay.pom.v.maom1 := sc.poc1 - sc.maoc1
    clay.pom.v.maom2 := sc.poc2 - sc.maoc1
    '

fit_model_SEM <- sem(model_SEM,
                     data = dat %>%
                       mutate(ET = ET/100,
                              prcp_mean = prcp_mean/100,
                              NPP = NPP*10,
                              al_fe = al_fe*10),
                     std.ov = FALSE, #do not standardize all variables before analysis
                     fixed.x = FALSE,
                     group = "AI_factor")

fit_model_SEM_equal <- sem(model_SEM,
                           data = dat %>%
                             mutate(ET = ET/100,
                                    prcp_mean = prcp_mean/100,
                                    NPP = NPP*10,
                                    al_fe = al_fe*10),
                           std.ov = FALSE, #do not standardize all variables before analysis
                           fixed.x = FALSE,
                           group = "AI_factor",
                           group.equal = c("intercepts", "regressions"))


anova(fit_model_SEM, fit_model_SEM_equal) #Multi-group model with some constraints is better than the constrained/equal model
summary(fit_model_SEM, fit.measures = TRUE, standardized = T, rsq = T)


#Figure 2A created outside of R
cat("Fig 2. N = ", summary(fit_model_SEM)$data$nobs %>% sum(), 
    "\n(Omitted", summary(fit_model_SEM)$data$norig %>% sum() - summary(fit_model_SEM)$data$nobs %>% sum(), "due to missing data)")

## Fig. 3B ----
sem_coeffs <- read.csv("Indirect paths.csv")

library(ggpattern)
plot_sem_coeffs_Arid <- ggplot(sem_coeffs %>%
                                 filter(Climate == "Arid") %>%
                                 filter(!(Med_variable %in% c("Indirect via pH and NPP", "Indirect via just NPP"))) %>%
                                 mutate(`SOC fraction` = factor(Response, levels=c('POC', 'MAOC')),
                                        pred_var_f = factor(Variable, levels=c('Cultivation', 'Depth', 'Caex + Mgex', 'Alo + [1/2]Feo', 'Silt + Clay', 'ET0', 'MAP', 'NPP', 'pH'))),
                               aes(y = Value, x = Med_variable, 
                                   fill = pred_var_f)) +
  geom_bar(position = "stack", stat = "identity") +
  facet_grid(cols = vars(`SOC fraction`)) +
  labs(y = "Standardized effect size", x = "") +
  scale_x_discrete(labels = c("Direct", expression(atop("Indirect", "via NPP")), expression(atop("Indirect", "via pH")))) +
  geom_hline(yintercept = 0) +
  theme_classic() +
  theme(legend.position="bottom") +
  scale_fill_manual(name = "Variable", values = c("#022d25", "#f4f4f5", "#ead8c3", "#cc8647", "#E1C696", "#9cd9e1", "#00A6AC", "#106d5c", "#C6C6C5"))

plot_sem_coeffs_Humid <- ggplot(sem_coeffs %>%
                                  filter(Climate == "Humid") %>%
                                  filter(!(Med_variable %in% c("Indirect via pH and NPP", "Indirect via just NPP"))) %>%
                                  mutate(`SOC fraction` = factor(Response, levels=c('POC', 'MAOC')),
                                         pred_var_f = factor(Variable, levels=c('Cultivation', 'Depth', 'Caex + Mgex', 'Alo + [1/2]Feo', 'Silt + Clay', 'ET0', 'MAP', 'NPP', 'pH'))),
                                aes(y = Value, x = Med_variable, 
                                    fill = pred_var_f)) +
  geom_bar(position = "stack", stat = "identity") +
  facet_grid(cols = vars(`SOC fraction`)) +
  labs(y = "Standardized effect size", x = "") +
  scale_x_discrete(labels = c("Direct", expression(atop("Indirect", "via NPP")), expression(atop("Indirect", "via pH")))) +
  geom_hline(yintercept = 0) +
  theme_classic() +
  theme(legend.position="bottom") +
  scale_fill_manual(name = "Variable", values = c("#022d25", "#f4f4f5", "#ead8c3", "#cc8647", "#E1C696", "#9cd9e1", "#00A6AC", "#106d5c", "#C6C6C5"))

#library(patchwork)
p3 <- (plot_sem_coeffs_Arid | plot_sem_coeffs_Humid) /
  guide_area() +
  plot_layout(guides = 'collect', heights = c(4.5, 1))
pdf('Fig_3B.pdf', 
    width = 10.1, 
    height= 5, 
    useDingbats=FALSE)
p3
dev.off()

# 4) Random Forest ----
library(randomForest)
library(caret)

# Create datasets for each category
vars_rf <- c("hzn_bot", "horizon_ap", 
             "ph_h2o", "ca_mg", "al_fe",
             "claysilt", "ET", "prcp_mean", "NPP") 
vars_C <- c("MAOM_mgpg_OC", "POM_mgpg_OC")

dat_rf_arid <- dat %>%
  filter(AI_factor == "Arid") %>%
  dplyr::select(labsampnum, all_of(c(vars_C, vars_rf))) %>% # this is dropping ~50 soils, consider using fewer columns to drop less
  mutate_if(is.character, as.factor) %>%
  drop_na()

dat_rf_humid <- dat %>%
  filter(AI_factor == "Humid") %>%
  dplyr::select(labsampnum, all_of(c(vars_C, vars_rf))) %>% # this is dropping ~50 soils, consider using fewer columns to drop less
  mutate_if(is.character, as.factor) %>%
  drop_na()
  
  
### MAOC Arid----
## Split the data so that we use 75% of it for training
## Don't need to cross-validate for random forest models, see https://stackoverflow.com/questions/31637259/random-forest-crossvalidation-in-r
set.seed(123)
train_index_arid <- createDataPartition(y = dat_rf_arid$MAOM_mgpg_OC,
                                        p = 0.75, list = FALSE) #Using 75% of data for training
training_set_MAOMC_arid <- dat_rf_arid[train_index_arid, ]  %>%
  dplyr::select(-any_of(vars_C[! vars_C %in% "MAOM_mgpg_OC"]), -labsampnum)
testing_set_MAOMC_arid <- as.data.frame(dat_rf_arid[-train_index_arid, ]) %>%
  dplyr::select(-any_of(vars_C[! vars_C %in% "MAOM_mgpg_OC"]), -labsampnum)

set.seed(123)
RF_mod_MAOC_arid = randomForest(MAOM_mgpg_OC ~ ., 
                                xtest = testing_set_MAOMC_arid[,-1],
                                ytest = testing_set_MAOMC_arid[,1],
                                data = training_set_MAOMC_arid,
                                ntree = 500,
                                localImp = TRUE,
                                type = "regression")
RF_mod_MAOC_arid # test set R2 = 52.15
plot(RF_mod_MAOC_arid, main = "Arid MAOM RF Model")
varImpPlot(RF_mod_MAOC_arid, main = "Arid MAOM RF Model, 53.19% var expl")

#### Look at predictions vs. actual values for testing set----
dat_MAOMC_pred_arid <- RF_mod_MAOC_arid$predicted %>%
  as.data.frame() %>%
  rename("MAOM_mgpg_OC_pred" = ".") %>%
  bind_cols(dat_rf_arid[train_index_arid, ]  %>% # training set predictions
              dplyr::select("MAOM_mgpg_OC_meas" = MAOM_mgpg_OC, labsampnum)) %>%
  mutate(set = "training") %>%
  bind_rows(RF_mod_MAOC_arid$test$predicted %>% #testing set predictions
              as.data.frame() %>%
              rename("MAOM_mgpg_OC_pred" = ".") %>%
              bind_cols(dat_rf_arid[-train_index_arid, ]  %>%
                          dplyr::select("MAOM_mgpg_OC_meas" = MAOM_mgpg_OC, labsampnum)) %>%
              mutate(set = "testing"))

ggplot(aes(x = MAOM_mgpg_OC_meas,
           y = MAOM_mgpg_OC_pred), data = dat_MAOMC_pred_arid) +
  geom_point(aes(fill = set), 
             alpha = 0.7, 
             size = 2, 
             shape = 21) +
  scale_fill_manual(values = c("blue", "red")) +
  scale_y_continuous(limits = c(0, 110)) +
  geom_abline(intercept = 0, slope = 1) +
  theme_classic() +
  labs(title = "Arid",
       x= expression(paste("Observed MAOM-C (mg g"^"-1", " whole soil)" )),
       y = expression(paste("Predicted MAOM-C (mg g"^"-1", " whole soil)" )))

### MAOC Humid ----
set.seed(123)
train_index_humid <- createDataPartition(y = dat_rf_humid$MAOM_mgpg_OC,
                                         p = 0.75, list = FALSE) #Using 75% of data for training
training_set_MAOMC_humid <- dat_rf_humid[train_index_humid, ]  %>%
  dplyr::select(-any_of(vars_C[! vars_C %in% "MAOM_mgpg_OC"]), -labsampnum)
testing_set_MAOMC_humid <- as.data.frame(dat_rf_humid[-train_index_humid, ]) %>%
  dplyr::select(-any_of(vars_C[! vars_C %in% "MAOM_mgpg_OC"]), -labsampnum)

set.seed(123)
RF_mod_MAOC_humid = randomForest(MAOM_mgpg_OC ~ ., 
                                 xtest = testing_set_MAOMC_humid[,-1],
                                 ytest = testing_set_MAOMC_humid[,1],
                                 data = training_set_MAOMC_humid,
                                 ntree = 500,
                                 localImp = TRUE,
                                 type = "regression")
RF_mod_MAOC_humid 
plot(RF_mod_MAOC_humid, main = "humid MAOM RF Model")
varImpPlot(RF_mod_MAOC_humid, main = "humid MAOM RF Model, 47.62% var expl")


#### Look at predictions vs. values for testing set ----
dat_MAOMC_pred_humid <- RF_mod_MAOC_humid$predicted %>%
  as.data.frame() %>%
  rename("MAOM_mgpg_OC_pred" = ".") %>%
  bind_cols(dat_rf_humid[train_index_humid, ]  %>% # training set predictions
              dplyr::select("MAOM_mgpg_OC_meas" = MAOM_mgpg_OC, labsampnum)) %>%
  mutate(set = "training") %>%
  bind_rows(RF_mod_MAOC_humid$test$predicted %>% #testing set predictions
              as.data.frame() %>%
              rename("MAOM_mgpg_OC_pred" = ".") %>%
              bind_cols(dat_rf_humid[-train_index_humid, ]  %>%
                          dplyr::select("MAOM_mgpg_OC_meas" = MAOM_mgpg_OC, labsampnum)) %>%
              mutate(set = "testing"))

ggplot(aes(x = MAOM_mgpg_OC_meas,
           y = MAOM_mgpg_OC_pred), data = dat_MAOMC_pred_humid) +
  geom_point(aes(fill = set), 
             alpha = 0.7, 
             size = 2, 
             shape = 21) +
  scale_fill_manual(values = c("blue", "red")) +
  scale_y_continuous(limits = c(0, 110)) +
  geom_abline(intercept = 0, slope = 1) +
  theme_classic() +
  labs(title = "Humid",
    x= expression(paste("Observed MAOM-C (mg g"^"-1", " whole soil)" )),
       y = expression(paste("Predicted MAOM-C (mg g"^"-1", " whole soil)" )))

### POC ----
### Arid
set.seed(123)
train_index_POC_arid <- createDataPartition(y = dat_rf_arid$POM_mgpg_OC,
                                            p = 0.75, list = FALSE) #Using 75% of data for training
training_set_POC_arid <- dat_rf_arid[train_index_POC_arid, ]  %>%
  dplyr::select(-any_of(vars_C[! vars_C %in% "POM_mgpg_OC"]), -labsampnum)
testing_set_POC_arid <- as.data.frame(dat_rf_arid[-train_index_POC_arid, ]) %>%
  dplyr::select(-any_of(vars_C[! vars_C %in% "POM_mgpg_OC"]), -labsampnum)

set.seed(123)
RF_mod_POC_arid = randomForest(POM_mgpg_OC ~ ., 
                               xtest = testing_set_POC_arid[,-1],
                               ytest = testing_set_POC_arid[,1],
                               data = training_set_POC_arid,
                               ntree = 500,
                               localImp = TRUE,
                               type = "regression")
RF_mod_POC_arid # test set R2 is good - 77%
plot(RF_mod_POC_arid, main = "Arid POM RF Model")
varImpPlot(RF_mod_POC_arid, main = "Arid POM RF Model, 33.79% var expl")

### Humid
set.seed(123)
train_index_POC_humid <- createDataPartition(y = dat_rf_humid$POM_mgpg_OC,
                                             p = 0.75, list = FALSE) #Using 75% of data for training
training_set_POC_humid <- dat_rf_humid[train_index_POC_humid, ]  %>%
  dplyr::select(-any_of(vars_C[! vars_C %in% "POM_mgpg_OC"]), -labsampnum)
testing_set_POC_humid <- as.data.frame(dat_rf_humid[-train_index_POC_humid, ]) %>%
  dplyr::select(-any_of(vars_C[! vars_C %in% "POM_mgpg_OC"]), -labsampnum)

set.seed(123)
RF_mod_POC_humid = randomForest(POM_mgpg_OC ~ ., 
                                xtest = testing_set_POC_humid[,-1],
                                ytest = testing_set_POC_humid[,1],
                                data = training_set_POC_humid,
                                ntree = 500,
                                localImp = TRUE,
                                type = "regression")
RF_mod_POC_humid 
plot(RF_mod_POC_humid, main = "humid POM RF Model")
varImpPlot(RF_mod_POC_humid, main = "humid POM RF Model, 30.2% var expl")

### Look at predictions vs. actual values for testing set----
# performs poorly at higher values.
#arid
dat_POMC_arid_pred <- RF_mod_POC_arid$predicted %>%
  as.data.frame() %>%
  rename("POM_mgpg_OC_pred" = ".") %>%
  bind_cols(dat_rf_arid[train_index_POC_arid, ]  %>% # training set predictions
              dplyr::select("POM_mgpg_OC_meas" = POM_mgpg_OC, labsampnum)) %>%
  mutate(set = "training") %>%
  bind_rows(RF_mod_POC_arid$test$predicted %>% #testing set predictions
              as.data.frame() %>%
              rename("POM_mgpg_OC_pred" = ".") %>%
              bind_cols(dat_rf_arid[-train_index_POC_arid, ]  %>%
                          dplyr::select("POM_mgpg_OC_meas" = POM_mgpg_OC, labsampnum)) %>%
              mutate(set = "testing"))

ggplot(aes(x = POM_mgpg_OC_meas,
           y = POM_mgpg_OC_pred), data = dat_POMC_arid_pred) +
  geom_point(aes(fill = set), 
             alpha = 0.7, 
             size = 2, 
             shape = 21) +
  scale_fill_manual(values = c("blue", "red")) +
  scale_y_continuous(limits = c(0, 110)) +
  geom_abline(intercept = 0, slope = 1) +
  theme_classic() +
  labs(title = "Arid",
    x= expression(paste("Observed POM-C (mg g"^"-1", " whole soil)" )),
       y = expression(paste("Predicted POM-C (mg g"^"-1", " whole soil)" )))

#humid
dat_POMC_humid_pred <- RF_mod_POC_humid$predicted %>%
  as.data.frame() %>%
  rename("POM_mgpg_OC_pred" = ".") %>%
  bind_cols(dat_rf_humid[train_index_POC_humid, ]  %>% # training set predictions
              dplyr::select("POM_mgpg_OC_meas" = POM_mgpg_OC, labsampnum)) %>%
  mutate(set = "training") %>%
  bind_rows(RF_mod_POC_humid$test$predicted %>% #testing set predictions
              as.data.frame() %>%
              rename("POM_mgpg_OC_pred" = ".") %>%
              bind_cols(dat_rf_humid[-train_index_POC_humid, ]  %>%
                          dplyr::select("POM_mgpg_OC_meas" = POM_mgpg_OC, labsampnum)) %>%
              mutate(set = "testing"))

ggplot(aes(x = POM_mgpg_OC_meas,
           y = POM_mgpg_OC_pred), data = dat_POMC_humid_pred) +
  geom_point(aes(fill = set), 
             alpha = 0.7, 
             size = 2, 
             shape = 21) +
  scale_fill_manual(values = c("blue", "red")) +
  scale_y_continuous(limits = c(0, 110)) +
  geom_abline(intercept = 0, slope = 1) +
  theme_classic() +
  labs(title = "Humid",
       x= expression(paste("Observed POM-C (mg g"^"-1", " whole soil)" )),
       y = expression(paste("Predicted POM-C (mg g"^"-1", " whole soil)" )))

### Fig. 4  ----
#### POC arid
imp_poc_arid <- varImpPlot(RF_mod_POC_arid) 
# create the data.frame for the plot 
p_imp_poc_arid <- imp_poc_arid %>% # save the varImp object
  as.data.frame() %>%
  rownames_to_column("varnames") %>%
  mutate(varnames = fct_recode(varnames, Depth = "hzn_bot", Cultivation = "horizon_ap",
                               pH = "ph_h2o", 'Caex + Mgex' = "ca_mg", 'Alo + [1/2]Feo' = "al_fe",
                               'Silt + Clay' = "claysilt", ET0 = "ET", MAP = "prcp_mean", NPP = "NPP"))

p4_imp_poc_arid <- ggplot(p_imp_poc_arid, aes(y = reorder(varnames, `%IncMSE`), weight = `%IncMSE`, fill=as.factor(varnames))) + 
  geom_bar() +
  scale_fill_manual(name="Variable", values = c("#cc8647", "#ead8c3","#E1C696", "#9cd9e1","#022d25", "#f4f4f5","#106d5c","#C6C6C5","#00A6AC")) +
  xlab("%IncMSE") +
  ylab("Variable") +
  theme_classic()

#### POC humid
imp_poc_humid <- varImpPlot(RF_mod_POC_humid) 
# create the data.frame for the plot 
p_imp_poc_humid <- imp_poc_humid %>% # save the varImp object
  as.data.frame() %>%
  rownames_to_column("varnames") %>%
  mutate(varnames = fct_recode(varnames, Depth = "hzn_bot", Cultivation = "horizon_ap",
                               pH = "ph_h2o", 'Caex + Mgex' = "ca_mg", 'Alo + [1/2]Feo' = "al_fe",
                               'Silt + Clay' = "claysilt", ET0 = "ET", MAP = "prcp_mean", NPP = "NPP"))

p4_imp_poc_humid <- ggplot(p_imp_poc_humid, aes(y = reorder(varnames, `%IncMSE`), weight = `%IncMSE`, fill=as.factor(varnames))) + 
  geom_bar() +
  scale_fill_manual(name="Variable", values = c("#cc8647", "#ead8c3","#E1C696", "#9cd9e1","#022d25", "#f4f4f5","#106d5c","#C6C6C5","#00A6AC")) +
  xlab("%IncMSE") +
  ylab("Variable") +
  theme_classic()


#### MAOC arid
imp_maoc_arid <- varImpPlot(RF_mod_MAOC_arid) 
# create the data.frame for the plot 
p_imp_maoc_arid <- imp_maoc_arid %>% # save the varImp object
  as.data.frame() %>%
  rownames_to_column("varnames") %>%
  mutate(varnames = fct_recode(varnames, Depth = "hzn_bot", Cultivation = "horizon_ap",
                               pH = "ph_h2o", 'Caex + Mgex' = "ca_mg", 'Alo + [1/2]Feo' = "al_fe",
                               'Silt + Clay' = "claysilt", ET0 = "ET", MAP = "prcp_mean", NPP = "NPP"))

p4_imp_maoc_arid <- ggplot(p_imp_maoc_arid, aes(y = reorder(varnames, `%IncMSE`), weight = `%IncMSE`, fill=as.factor(varnames))) + 
  geom_bar() +
  scale_fill_manual(name="Variable", values = c("#cc8647", "#ead8c3","#E1C696", "#9cd9e1","#022d25", "#f4f4f5","#106d5c","#C6C6C5","#00A6AC")) +
  xlab("%IncMSE") +
  ylab("Variable") +
  theme_classic()

#### MAOC arid
imp_maoc_humid <- varImpPlot(RF_mod_MAOC_humid) 
# create the data.frame for the plot 
p_imp_maoc_humid <- imp_maoc_humid %>% # save the varImp object
  as.data.frame() %>%
  rownames_to_column("varnames") %>%
  mutate(varnames = fct_recode(varnames, Depth = "hzn_bot", Cultivation = "horizon_ap",
                               pH = "ph_h2o", 'Caex + Mgex' = "ca_mg", 'Alo + [1/2]Feo' = "al_fe",
                               'Silt + Clay' = "claysilt", ET0 = "ET", MAP = "prcp_mean", NPP = "NPP"))

p4_imp_maoc_humid <- ggplot(p_imp_maoc_humid, aes(y = reorder(varnames, `%IncMSE`), weight = `%IncMSE`, fill=as.factor(varnames))) + 
  geom_bar() +
  scale_fill_manual(name="Variable", values = c("#cc8647", "#ead8c3","#E1C696", "#9cd9e1","#022d25", "#f4f4f5","#106d5c","#C6C6C5","#00A6AC")) +
  xlab("%IncMSE") +
  ylab("Variable") +
  theme_classic()

### Make plot 
p4 <- p4_imp_poc_arid + p4_imp_poc_humid + p4_imp_maoc_arid + p4_imp_maoc_humid +
  plot_layout(guides = 'collect') +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 16))

pdf('Fig_4.pdf', 
    width = 10, 
    height= 6, 
    useDingbats=FALSE)
p4
dev.off()
# 5) Saturation ----
library(quantreg)
library(visreg)

rqfit_1 <- rq(MAOM_mgpg_OC ~ claysilt:AI_factor - 1, 
              tau = 0.95, data = dat)
summary(rqfit_1)
visreg(rqfit_1, "claysilt", by = "AI_factor", overlay=TRUE)


quant_reg_lines_clim <- data.frame(AI_factor = factor(c("Arid", "Humid"), levels = c("Arid", "Humid")),
                                   intercept = c(0,0),
                                   slope = summary(rqfit_1)$coefficients[1:2],
                                   se = summary(rqfit_1)$coefficients[3:4]) %>%
  mutate(ci_lower = slope - se*1.96,
         ci_upper = slope + se*1.96)

## Fig. 5 ----
Fig_5 <-  
  ggplot(dat,
         aes(x = claysilt, y = MAOM_mgpg_OC)) +
  geom_point(aes(color = AI_factor), alpha = 0.5, shape=16, size = 1.2) +
  geom_point(aes(x = claysilt, y = MAOM_mgpg_OC),
             dat %>% filter(WS_mgpg_OC < 100 & AI_factor == "Arid"),
             alpha = 0.5, color = "#eeb868", shape = 16, size = 1.2) +
  geom_abline(slope = seq(quant_reg_lines_clim$ci_lower[1],
                          quant_reg_lines_clim$ci_upper[1], 0.001),
              intercept = 0,
              color = "grey80", alpha = .15) +
  geom_abline(slope = seq(quant_reg_lines_clim$ci_lower[2],
                          quant_reg_lines_clim$ci_upper[2], 0.001),
              intercept = 0,
              color = "grey80", alpha = 0.1) +
  xlab(expression(paste("Silt and clay content (%)"))) +
  ylab(expression(paste("MAOM C (g C kg"^-1, "soil)")))+
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(limits = c(0, 90), expand = c(0, 0)) +
  geom_abline(aes(slope = slope, intercept = intercept, 
                  color = AI_factor, linetype = AI_factor),
              data = quant_reg_lines_clim, linewidth = .7) +
  theme_classic() +
  scale_fill_manual(name = "Aridity Group",
                    values = c( "#eeb868", "#49beaa")) +
  scale_color_manual(name = "Aridity Group", values = c("#eeb868", "#49beaa")) +
  scale_linetype_manual(name = "Aridity Group", values = c(1, 2, 6)) +
  geom_abline(slope = 0.86, intercept = 0, color = "gray60",
               linewidth = .7, linetype = 1) +
  theme(legend.position = c(.2, .85))

cat("Fig 5. N = ", nrow(filter(dat, !is.na(claysilt))), "\n directly measured = ",
    nrow(filter(dat, !is.na(claysilt) & fraction_quant_method == "direct measurement")),
    "MIR predicted = ", nrow(filter(dat, !is.na(claysilt) & fraction_quant_method == "MIR prediction")))

pdf('Fig_5.pdf', 
    width = 4, 
    height= 4, 
    useDingbats=FALSE)
Fig_5
dev.off()


# Fig. S6 ----
library(ggpubr)    
Fig_S6A <- 
  ggplot(dat, 
         aes(x = log(SOC_NPP_Ratio), y = fMAOM)) +
  geom_point(aes(shape = land_cover_simple), alpha = .4, color = "#6a9abb") +
  scale_shape_manual(name = "Land cover", values = c(4, 16, 2, 8, 3)) +
  guides(shape = "none") +
  geom_smooth(method = "lm", linewidth = .7, color = "#064672") +
  stat_cor(aes(label = paste(after_stat(rr.label), after_stat(p.label), 
                             sep = "~`,`~")),
           label.x.npc = 0.05, label.y.npc = 0.1, show.legend = FALSE) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1.06)) +
  labs(color = "Land Cover", x = expression(paste("SOC:NPP ratio")),
       y = expression(f[MAOC])) +
  theme_classic() 

Fig_S6B <- 
  ggplot(aes(x = AI_factor), dat = dat) +
  geom_point(aes(y = fMAOM, shape = land_cover_simple), color = "#6a9abb", alpha = .4, position = "jitter",
             size = 2) +
  geom_boxplot(aes(y = fMAOM), color = "#064672", alpha = 0) +
  scale_shape_manual(name = "Land cover", values = c(4, 16, 2, 8, 3)) +
  #guides(shape = "none") +
  labs(x = "Aridity Index", y =  
         expression(f[MAOC])) +
  theme_classic() +
  theme(axis.title.y = element_text(color = "#064672"),
        legend.position = "none") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.08)) +
  annotate("text", x = c(1, 2), y=c(rep(1.05, 2)), label = c("a", "b"))
  

#fMAOC vs. SOC Storage Efficiency, by land cover
f_labels <- data.frame(land_cover_simple = c("cropland", "forest", "grassland", "other", "shrubland"),
                       label = c("*", "*", "*", "", "")) #not significant for "other" and "shrubland"

Fig_S6C <- ggplot(data = dat %>%
                                      mutate(land_cover_simple = factor(land_cover_simple),
                                             horizon_ap = forcats::fct_relevel(horizon_ap, "1", "0"),
                                             SOC_NPP_Ratio_factor = forcats::fct_relevel(SOC_NPP_Ratio_factor, "low", "high"),
                                             SOC_NPP_Ratio_factor = forcats::fct_recode(SOC_NPP_Ratio_factor, "Low" = "low", "High" = "high")),
                                    aes(x = SOC_NPP_Ratio_factor, y = fMAOM)) +
  geom_point(aes(shape = land_cover_simple, color = SOC_NPP_Ratio_factor), alpha = 0.4,
             position = position_jitterdodge()) +
  #No support in general ANOVA for different groupings of land cover. (no land cover : storage_eff intxn)
  geom_boxplot(aes(color = SOC_NPP_Ratio_factor), alpha = 0) +
  scale_shape_manual(name = "Land cover", values = c(4, 16, 2, 8, 3)) +
  scale_color_manual(name = "SOC:NPP", values = c("#064672", "#6a9abb")) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1.06)) +
  facet_wrap(~ land_cover_simple, strip.position = "bottom", nrow = 1) +
  labs(shape = "Land Cover", x = NULL, y = expression(f[MAOC]))+
  theme_classic() +
  theme(strip.placement = "outside") +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x = element_blank()) +
  guides(shape = guide_legend(override.aes = list(color = "black", size = 2, alpha = .8))) +
  geom_text(x = 1.5, y = 1.01, aes(label = label), data = f_labels, size = 6)

# independent 2-group Mann-Whitney U Test, by land cover
wilcox.test(fMAOM ~ SOC_NPP_Ratio_factor, data = dat %>%
              filter(land_cover_simple == "cropland")) #W = 14806, p-value = 0.00144
wilcox.test(fMAOM ~ SOC_NPP_Ratio_factor, data = dat %>%
              filter(land_cover_simple == "forest")) #W = 15930, p-value = 4.393e-10
wilcox.test(fMAOM ~ SOC_NPP_Ratio_factor, data = dat %>%
              filter(land_cover_simple == "grassland")) #W = 6000, p-value = 0.03148
wilcox.test(fMAOM ~ SOC_NPP_Ratio_factor, data = dat %>%
              filter(land_cover_simple == "other")) #W = 15, p-value = 0.152
wilcox.test(fMAOM ~ SOC_NPP_Ratio_factor, data = dat %>%
              filter(land_cover_simple == "shrubland")) #W = 48, p-value = 0.3364

# Comparing fMAOM between arid and humid
kruskal.test(fMAOM ~ AI_factor, data = dat) # p = 0.018


pS6 <- (Fig_S6A + Fig_S6B) /
  Fig_S6C + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 16))

pdf('Fig_S6.pdf', 
    width = 9, 
    height= 10, 
    useDingbats=FALSE)
pS6
dev.off()
