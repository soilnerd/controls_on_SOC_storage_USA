### Code for Lavallee et al. 2024 ###
# 0) Read data ----
library(tidyverse)
setwd("/Users/jocelynlavallee/Dropbox/3\ Other\ Projects/NRCS\ Soil\ Inventory/Manuscript")

load("Lavallee_et_al_2024_data.RData")
load("Lavallee_et_al_2024_supp_data.RData")


# 1) SOC and NPP vs. Aritidy Index ----
Fig_1A <-  
  ggplot(data = dat,
         aes(x = aridity_index_v3)) +
  geom_point(aes(y = WS_mgpg_OC_calculated*0.01), color = "#E1BE6A", 
             alpha = .7) +
  geom_point(data = dat,
             aes(y = NPP_MODIS_mean_2001_2022_kgC_m2), alpha = .7, col = "#7ec4a6") +
  xlab(expression(paste("Aridity Index"))) +
  geom_smooth(aes(y = NPP_MODIS_mean_2001_2022_kgC_m2), color = "#299d6c", method = "gam") +
  geom_smooth(aes(y = WS_mgpg_OC_calculated*0.01), color = "#b49854", method = "gam") +
  scale_y_continuous(name = expression(paste("NPP (kg C m"^-2, ")")),
                     sec.axis = sec_axis(~./.01, name = expression(paste("SOC (g C kg"^-1, "soil)")))) +
  scale_x_continuous(limits = c(0,4)) + # Omits 4 points > 4. Note the SEM analysis omitted points > 2.65
  theme_classic() +
  theme(axis.title.y = element_text(color = "#299d6c"),
        axis.title.y.right = element_text(color = "#b49854")) +
  geom_vline(xintercept = 0.5, color = "#7a7a7a") +
  geom_vline(xintercept = 0.9, color = "#7a7a7a")



Fig_1B <- 
  ggplot(dat, aes(x = aridity_index_v3)) +
  geom_point(aes(y = SOC_NPP_Ratio), color = "#6a9abb", alpha = .7) +
  geom_smooth(aes(y = SOC_NPP_Ratio), color = "#064672", method = "gam") +
  labs(x = "Aridity Index", y =  
         expression(paste("SOC:NPP ratio"))) +
  scale_x_continuous(limits = c(0,4)) + # Omits 4 points > 4. Note the SEM analysis omitted points > 2.65
  theme_classic() +
  theme(axis.title.y = element_text(color = "#064672")) +
  geom_vline(xintercept = 0.5, color = "#7a7a7a")+
  geom_vline(xintercept = 0.9, color = "#7a7a7a")

kruskal.test(SOC_NPP_Ratio ~ aridity_index_factor, data = dat)
library(FSA)
dunnTest(SOC_NPP_Ratio ~ aridity_index_factor, data = dat,
         method = "holm")

cat("Fig 1A and B  N = ", nrow(filter(dat, aridity_index_v3 < 4)))

Fig_1C <- 
  ggplot(aes(x = aridity_index_factor), dat = dat) +
  geom_point(aes(y = log(SOC_NPP_Ratio)), color = "#6a9abb", alpha = .2, position = "jitter",
             size = 1) +
  geom_boxplot(aes(y = log(SOC_NPP_Ratio)), color = "#064672", alpha = 0) +
  labs(x = "Aridity Index", y =  
         expression(paste("ln(SOC:NPP ratio)"))) +
  theme_classic() +
  theme(axis.title.y = element_text(color = "#064672")) +
  scale_y_continuous(expand = c(0, 0), limits = c(1, 6.5)) + # Omits 16 points outside range for readability
  annotate("text", x = c(1, 2, 3), y=c(rep(6.3, 3)), label = c("a", "b", "c"))

cat("Fig 1C  N = ", nrow(filter(dat, log(SOC_NPP_Ratio) > 1))) # no points > 6.5


library(patchwork)
p1 <- Fig_1A + Fig_1B + inset_element(Fig_1C, 0.45, 0.45, 1, 1.05) +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 16))
pdf('Fig_1.pdf', 
    width = 12, 
    height= 5.5, 
    useDingbats=FALSE)
p1
dev.off()


# 2) Structural Equation Modeling ----
library(lavaan)
model_SEM <- ' 
 POM_mgpg_OC ~ 
      c("clay.poc1", "clay.poc1", "clay.poc1") * clay_total + 
      c("ph.poc1", "ph.poc1", "ph.poc2") * ph_h2o + 
      c("mci.poc1", "mci.poc2", "mci.poc2") * MCI + 
      c("et.poc1", "et.poc2", "et.poc1") * ET + 
      c("map.poc1", "map.poc2", "map.poc1") * prcp_mean +  
      c("npp.poc1", "npp.poc2", "npp.poc2") * NPP_MODIS_mean_2001_2022_kgC_m2 + 
      c("cult.poc1", "cult.poc1", "cult.poc2") * horizon_ap + 
      c("hzn.poc1", "hzn.poc1", "hzn.poc1") * hzn_bot 
    MAOM_mgpg_OC ~ 
      c("clay.maoc1", "clay.maoc1", "clay.maoc1") * clay_total + 
      c("ph.maoc1", "ph.maoc2", "ph.maoc2") * ph_h2o + 
      c("mci.maoc1", "mci.maoc2", "mci.maoc3") * MCI + 
      c("et.maoc1", "et.maoc2", "et.maoc2") * ET +
      c("map.maoc1", "map.maoc2", "map.maoc2") * prcp_mean + 
      c("npp.maoc1", "npp.maoc2", "npp.maoc2") * NPP_MODIS_mean_2001_2022_kgC_m2 + 
      c("cult.maoc1", "cult.maoc1", "cult.maoc1") *  horizon_ap + 
      c("hzn.maoc1", "hzn.maoc1", "hzn.maoc1") * hzn_bot 
    ph_h2o ~ 
      c("et.ph1", "et.ph2", "et.ph2") * ET + 
      c("map.ph1", "map.ph2", "map.ph3") * prcp_mean + 
      c("mci.ph1", "mci.ph2", "mci.ph3") * MCI + 
      c("cult.ph1", "cult.ph2", "cult.ph2") * horizon_ap + 
      c("hzn.ph1", "hzn.ph1", "hzn.ph1") * hzn_bot 
    NPP_MODIS_mean_2001_2022_kgC_m2 ~ 
      c("clay.npp1", "clay.npp2", "clay.npp2") * clay_total +  
      c("ph.npp1", "ph.npp1", "ph.npp1") * ph_h2o + 
      c("et.npp1", "et.npp2", "et.npp3") * ET + 
      c("map.npp1", "map.npp2", "map.npp3") * prcp_mean + 
      c("mci.npp1", "mci.npp2", "mci.npp1") * MCI + 
      c("cult.npp1", "cult.npp2", "cult.npp2") * horizon_ap + 
      c("hzn.npp1", "hzn.npp1", "hzn.npp2") * hzn_bot 

 # direct paths
              direct.clay.poc1 := clay.poc1
                  direct.clay.poc2 := clay.poc1
                  direct.clay.poc3 := clay.poc1
              direct.ph.poc1 := ph.poc1
                  direct.ph.poc2 := ph.poc1
                      direct.ph.poc3 := ph.poc2
              direct.mci.poc1 := mci.poc1
                      direct.mci.poc2 := mci.poc2
                      direct.mci.poc3 := mci.poc2
              direct.et.poc1 := et.poc1
                      direct.et.poc2 := et.poc2
                  direct.et.poc3 := et.poc1
              direct.map.poc1 := map.poc1
                      direct.map.poc2 := map.poc2
                  direct.map.poc3 := map.poc1
              direct.npp.poc1 := npp.poc1 
                      direct.npp.poc2 := 0
                      direct.npp.poc3 := 0
              direct.cult.poc1 := cult.poc1
                  direct.cult.poc2 := cult.poc1
                      direct.cult.poc3 := 0
              direct.hzn.poc1 := hzn.poc1
                  direct.hzn.poc2 := hzn.poc1
                  direct.hzn.poc3 := hzn.poc1
                  
              direct.clay.maoc1 := clay.maoc1
                   direct.clay.maoc2 := clay.maoc1
                   direct.clay.maoc3 := clay.maoc1
              direct.ph.maoc1 := ph.maoc1
                            direct.ph.maoc2 := ph.maoc2
                            direct.ph.maoc3 := ph.maoc2
              direct.mci.maoc1 := mci.maoc1
                 direct.mci.maoc2 := mci.maoc2
                   direct.mci.maoc3 := mci.maoc3
              direct.et.maoc1 := et.maoc1
                            direct.et.maoc2 := et.maoc2
                            direct.et.maoc3 := et.maoc2
              direct.map.maoc1 := map.maoc1
                            direct.map.maoc2 := 0
                            direct.map.maoc3 := 0
              direct.npp.maoc1 := npp.maoc1
                            direct.npp.maoc2 := 0
                            direct.npp.maoc3 := 0
              direct.cult.maoc1 := cult.maoc1
                   direct.cult.maoc2 := cult.maoc1
                   direct.cult.maoc3 := cult.maoc1
              direct.hzn.maoc1 := hzn.maoc1
                   direct.hzn.maoc2 := hzn.maoc1
                   direct.hzn.maoc3 := hzn.maoc1
                            
 #  indirect effects on SOC
           # Through pH (change/check pH paths)
             et.ph.poc1 := et.ph1 * direct.ph.poc1
                 et.ph.poc2 := et.ph2 * direct.ph.poc2
                 et.ph.poc3 := et.ph2 * direct.ph.poc3
             mci.ph.poc1 := mci.ph1 * direct.ph.poc1
                 mci.ph.poc2 := mci.ph2 * direct.ph.poc2
                   mci.ph.poc3 := mci.ph3 * direct.ph.poc3
             map.ph.poc1 := map.ph1 * direct.ph.poc1
                 map.ph.poc2 := map.ph2 * direct.ph.poc2 
                   map.ph.poc3 := map.ph3 * direct.ph.poc3 
             cult.ph.poc1 := 0 * direct.ph.poc1 
                 cult.ph.poc2 := cult.ph2 * direct.ph.poc2
                 cult.ph.poc3 := cult.ph2 * direct.ph.poc3
             hzn.ph.poc1 := hzn.ph1 * direct.ph.poc1
                 hzn.ph.poc2 := hzn.ph1 * direct.ph.poc2
                 hzn.ph.poc3 := hzn.ph1 * direct.ph.poc3
             
 
             et.ph.maoc1 := et.ph1 * direct.ph.maoc1
                 et.ph.maoc2 := et.ph2 * direct.ph.maoc2
                 et.ph.maoc3 := et.ph2 * direct.ph.maoc3
             mci.ph.maoc1 := mci.ph1 * direct.ph.maoc1
                 mci.ph.maoc2 := mci.ph2 * direct.ph.maoc2
                   mci.ph.maoc3 := mci.ph3 * direct.ph.maoc3
             map.ph.maoc1 := map.ph1 * direct.ph.maoc1
                 map.ph.maoc2 := map.ph2 * direct.ph.maoc2 
                   map.ph.maoc3 := map.ph3 * direct.ph.maoc3 
             cult.ph.maoc1 := 0 * direct.ph.maoc1 
                 cult.ph.maoc2 := cult.ph2 * direct.ph.maoc2
                 cult.ph.maoc3 := cult.ph2 * direct.ph.maoc3
             hzn.ph.maoc1 := hzn.ph1 * direct.ph.maoc1
                 hzn.ph.maoc2 := hzn.ph1 * direct.ph.maoc2
                 hzn.ph.maoc3 := hzn.ph1 * direct.ph.maoc3       
            
        # Through NPP  (change/check npp paths)
             clay.npp.poc1 := 0 * direct.npp.poc1
                    clay.npp.poc2 := clay.npp2 * direct.npp.poc2
                    clay.npp.poc3 := clay.npp2 * direct.npp.poc3
             ph.npp.poc1 := ph.npp1 * direct.npp.poc1
                 ph.npp.poc2 := ph.npp1 * direct.npp.poc2
                 ph.npp.poc3 := ph.npp1 * direct.npp.poc3
             et.npp.poc1 := et.npp1 * direct.npp.poc1
                    et.npp.poc2 := et.npp2 * direct.npp.poc2
                       et.npp.poc3 := et.npp3 * direct.npp.poc3
             mci.npp.poc1 := mci.npp1 * direct.npp.poc1
                    mci.npp.poc2 :=  0 * direct.npp.poc2
                 mci.npp.poc3 := mci.npp1 * direct.npp.poc3
             map.npp.poc1 := map.npp1 * direct.npp.poc1
                     map.npp.poc2 := 0 * direct.npp.poc2
                 map.npp.poc3 := map.npp1 * direct.npp.poc3
             cult.npp.poc1 := cult.npp1 * direct.npp.poc1
                  cult.npp.poc2 := cult.npp2 * direct.npp.poc2
                  cult.npp.poc3 := cult.npp2 * direct.npp.poc3
             hzn.npp.poc1 := hzn.npp1 * direct.npp.poc1
                hzn.npp.poc2 := hzn.npp1 * direct.npp.poc2
                     hzn.npp.poc3 := hzn.npp2 * direct.npp.poc3     
 
             clay.npp.maoc1 := 0 * direct.npp.maoc1
                    clay.npp.maoc2 := clay.npp2 * direct.npp.maoc2
                    clay.npp.maoc3 := clay.npp2 * direct.npp.maoc3
             ph.npp.maoc1 := ph.npp1 * direct.npp.maoc1
                 ph.npp.maoc2 := ph.npp1 * direct.npp.maoc2
                 ph.npp.maoc3 := ph.npp1 * direct.npp.maoc3
             et.npp.maoc1 := et.npp1 * direct.npp.maoc1
                    et.npp.maoc2 := et.npp2 * direct.npp.maoc2
                       et.npp.maoc3 := et.npp3 * direct.npp.maoc3
             mci.npp.maoc1 := mci.npp1 * direct.npp.maoc1
                    mci.npp.maoc2 :=  0 * direct.npp.maoc2
                 mci.npp.maoc3 := mci.npp1 * direct.npp.maoc3
             map.npp.maoc1 := map.npp1 * direct.npp.maoc1
                     map.npp.maoc2 := 0 * direct.npp.maoc2
                 map.npp.maoc3 := map.npp1 * direct.npp.maoc3
             cult.npp.maoc1 := cult.npp1 * direct.npp.maoc1
                  cult.npp.maoc2 := cult.npp2 * direct.npp.maoc2
                  cult.npp.maoc3 := cult.npp2 * direct.npp.maoc3
             hzn.npp.maoc1 := hzn.npp1 * direct.npp.maoc1
                hzn.npp.maoc2 := hzn.npp1 * direct.npp.maoc2
                     hzn.npp.maoc3 := hzn.npp2 * direct.npp.maoc3     
              
        # Through pH AND NPP 
             et.ph.npp.poc1 := et.ph1 * ph.npp1 * direct.npp.poc1
                 et.ph.npp.poc2 := et.ph1 * ph.npp1 * direct.npp.poc2
                 et.ph.npp.poc3 := et.ph2 * ph.npp1 * direct.npp.poc3
             mci.ph.npp.poc1 := mci.ph1 * ph.npp1 * direct.npp.poc1
                 mci.ph.npp.poc2 :=  mci.ph2 * ph.npp1 * direct.npp.poc2
                    mci.ph.npp.poc3 := mci.ph3 * ph.npp1 * direct.npp.poc3
             map.ph.npp.poc1 := map.ph1 * ph.npp1 * direct.npp.poc1
                 map.ph.npp.poc2 := map.ph2 * ph.npp1 * direct.npp.poc2
                   map.ph.npp.poc3 := map.ph3 * ph.npp1 * direct.npp.poc3
             cult.ph.npp.poc1 := 0 * ph.npp1 * direct.npp.poc1
                 cult.ph.npp.poc2 := cult.ph2 * ph.npp1 * direct.npp.poc2
              cult.ph.npp.poc3 := cult.ph2 * ph.npp1 * direct.npp.poc3
             hzn.ph.npp.poc1 := hzn.ph1 * ph.npp1 * direct.npp.poc1
                 hzn.ph.npp.poc2 := hzn.ph1 * ph.npp1 * direct.npp.poc2
              hzn.ph.npp.poc3 := hzn.ph1 * ph.npp1 * direct.npp.poc3
              
             et.ph.npp.maoc1 := et.ph1 * ph.npp1 * direct.npp.maoc1
                 et.ph.npp.maoc2 := et.ph1 * ph.npp1 * direct.npp.maoc2
                 et.ph.npp.maoc3 := et.ph2 * ph.npp1 * direct.npp.maoc3
             mci.ph.npp.maoc1 := mci.ph1 * ph.npp1 * direct.npp.maoc1
                 mci.ph.npp.maoc2 :=  mci.ph2 * ph.npp1 * direct.npp.maoc2
                   mci.ph.npp.maoc3 := mci.ph3 * ph.npp1 * direct.npp.maoc3
             map.ph.npp.maoc1 := map.ph1 * ph.npp1 * direct.npp.maoc1
                 map.ph.npp.maoc2 := map.ph2 * ph.npp1 * direct.npp.maoc2
                   map.ph.npp.maoc3 := map.ph3 * ph.npp1 * direct.npp.maoc3
             cult.ph.npp.maoc1 := 0 * ph.npp1 * direct.npp.maoc1
                 cult.ph.npp.maoc2 := cult.ph2 * ph.npp1 * direct.npp.maoc2
              cult.ph.npp.maoc3 := cult.ph2 * ph.npp1 * direct.npp.maoc3
             hzn.ph.npp.maoc1 := hzn.ph1 * ph.npp1 * direct.npp.maoc1
                 hzn.ph.npp.maoc2 := hzn.ph1 * ph.npp1 * direct.npp.maoc2
              hzn.ph.npp.maoc3 := hzn.ph1 * ph.npp1 * direct.npp.maoc3
             
  # total effects on SOC
             total.cult.poc1 := direct.cult.poc1 + cult.ph.poc1 + cult.npp.poc1 + cult.ph.npp.poc1
             total.cult.poc2 := direct.cult.poc2 + cult.ph.poc2 + cult.npp.poc2 + cult.ph.npp.poc2
             total.cult.poc3 := direct.cult.poc3 + cult.ph.poc3 + cult.npp.poc3 + cult.ph.npp.poc3
             total.et.poc1 := direct.et.poc1 + et.npp.poc1 + et.ph.poc1 + et.ph.npp.poc1
             total.et.poc2 := direct.et.poc2 + et.npp.poc2 + et.ph.poc2 + et.ph.npp.poc2
             total.et.poc3 := direct.et.poc3 + et.npp.poc3 + et.ph.poc3 + et.ph.npp.poc3
             total.mci.poc1 :=  direct.mci.poc1 + mci.ph.poc1 + mci.npp.poc1 + mci.ph.npp.poc1
             total.mci.poc2 :=  direct.mci.poc2 + mci.ph.poc2 + mci.npp.poc2 + mci.ph.npp.poc2
             total.mci.poc3 :=  direct.mci.poc3 + mci.ph.poc3 + mci.npp.poc3 + mci.ph.npp.poc3
             total.map.poc1 :=  direct.map.poc1 + map.ph.poc1 + map.npp.poc1 + map.ph.npp.poc1
             total.map.poc2 :=  direct.map.poc2 + map.ph.poc2 + map.npp.poc2 + map.ph.npp.poc2
             total.map.poc3 :=  direct.map.poc3 + map.ph.poc3 + map.npp.poc3 + map.ph.npp.poc3
             total.hzn.poc1 :=  direct.hzn.poc1 + hzn.ph.poc1 + hzn.npp.poc1 + hzn.ph.npp.poc1 
             total.hzn.poc2 :=  direct.hzn.poc2 + hzn.ph.poc2 + hzn.npp.poc2 + hzn.ph.npp.poc2
             total.hzn.poc3 :=  direct.hzn.poc3 + hzn.ph.poc3 + hzn.npp.poc3 + hzn.ph.npp.poc3 
             total.clay.poc1 := direct.clay.poc1 + clay.npp.poc1
             total.clay.poc2 := direct.clay.poc2 + clay.npp.poc2
             total.clay.poc3 := direct.clay.poc3 + clay.npp.poc3
             total.ph.poc1 := direct.ph.poc1 + ph.npp.poc1
             total.ph.poc2 := direct.ph.poc2 + ph.npp.poc2
             total.ph.poc3 := direct.ph.poc3 + ph.npp.poc3
    
             total.cult.maoc1 := direct.cult.maoc1 + cult.ph.maoc1 + cult.npp.maoc1 + cult.ph.npp.maoc1
             total.cult.maoc2 := direct.cult.maoc2 + cult.ph.maoc2 + cult.npp.maoc2 + cult.ph.npp.maoc2
             total.cult.maoc3 := direct.cult.maoc3 + cult.ph.maoc3 + cult.npp.maoc3 + cult.ph.npp.maoc3
             total.et.maoc1 := direct.et.maoc1 + et.npp.maoc1 + et.ph.maoc1 + et.ph.npp.maoc1
             total.et.maoc2 := direct.et.maoc2 + et.npp.maoc2 + et.ph.maoc2 + et.ph.npp.maoc2
             total.et.maoc3 := direct.et.maoc3 + et.npp.maoc3 + et.ph.maoc3 + et.ph.npp.maoc3
             total.mci.maoc1 :=  direct.mci.maoc1 + mci.ph.maoc1 + mci.npp.maoc1 + mci.ph.npp.maoc1
             total.mci.maoc2 :=  direct.mci.maoc2 + mci.ph.maoc2 + mci.npp.maoc2 + mci.ph.npp.maoc2
             total.mci.maoc3 :=  direct.mci.maoc3 + mci.ph.maoc3 + mci.npp.maoc3 + mci.ph.npp.maoc3
             total.map.maoc1 :=  direct.map.maoc1 + map.ph.maoc1 + map.npp.maoc1 + map.ph.npp.maoc1
             total.map.maoc2 :=  direct.map.maoc2 + map.ph.maoc2 + map.npp.maoc2 + map.ph.npp.maoc2
             total.map.maoc3 :=  direct.map.maoc3 + map.ph.maoc3 + map.npp.maoc3 + map.ph.npp.maoc3
             total.hzn.maoc1 :=  direct.hzn.maoc1 + hzn.ph.maoc1 + hzn.npp.maoc1 + hzn.ph.npp.maoc1 
             total.hzn.maoc2 :=  direct.hzn.maoc2 + hzn.ph.maoc2 + hzn.npp.maoc2 + hzn.ph.npp.maoc2
             total.hzn.maoc3 :=  direct.hzn.maoc3 + hzn.ph.maoc3 + hzn.npp.maoc3 + hzn.ph.npp.maoc3 
             total.clay.maoc1 := direct.clay.maoc1 + clay.npp.maoc1
             total.clay.maoc2 := direct.clay.maoc2 + clay.npp.maoc2
             total.clay.maoc3 := direct.clay.maoc3 + clay.npp.maoc3
             total.ph.maoc1 := direct.ph.maoc1 + ph.npp.maoc1
             total.ph.maoc2 := direct.ph.maoc2 + ph.npp.maoc2
             total.ph.maoc3 := direct.ph.maoc3 + ph.npp.maoc3
             
    # Differences, paths of interest
    # see https://groups.google.com/g/lavaan/c/I2Svn1gfTTU/m/Ma5Oouz2AAAJ
    cult.pom.v.maom12 := direct.cult.poc1 - direct.cult.maoc1 
    cult.pom.v.maom3 := direct.cult.poc3 - direct.cult.maoc3 
    hzn.pom.v.maom123 := direct.hzn.poc1 - direct.hzn.maoc1 
    ph.pom.v.maom1 := direct.ph.poc1 - direct.ph.maoc1 
    ph.pom.v.maom2 := direct.ph.poc1 - direct.ph.maoc2     
    ph.pom.v.maom3 := direct.ph.poc2 - direct.ph.maoc2 
    mci.pom.v.maom1 := direct.mci.poc1 - direct.mci.maoc1 
    mci.pom.v.maom2 := direct.mci.poc2 - direct.mci.maoc2
    mci.pom.v.maom3 := direct.mci.poc2- direct.mci.maoc3
    et.pom.v.maom1 := direct.et.poc1 - direct.et.maoc1
    et.pom.v.maom2 := direct.et.poc2 - direct.et.maoc2
    et.pom.v.maom3 := direct.et.poc1 - direct.et.maoc2
    map.pom.v.maom1 := direct.map.poc1 - direct.map.maoc1
    map.pom.v.maom2 := direct.map.poc2 - direct.map.maoc2
    map.pom.v.maom3 := direct.map.poc1 - direct.map.maoc2    
    npp.pom.v.maom1 := direct.npp.poc1 - direct.npp.maoc1
    clay.pom.v.maom123 := direct.clay.poc1 - direct.clay.maoc1
    '


fit_model_SEM <- sem(model_SEM,
                     data = dat %>%
                       mutate(ET = ET/100,
                              prcp_mean = prcp_mean/100,
                              NPP_MODIS_mean_2001_2022_kgC_m2 = NPP_MODIS_mean_2001_2022_kgC_m2*10),
                     std.ov = FALSE, 
                     fixed.x = FALSE,
                     group = "aridity_index_factor")


fit_model_SEM_equal <- sem(model_SEM,
                           data = dat %>%
                             mutate(ET = ET/100,
                                    prcp_mean = prcp_mean/100,
                                    NPP_MODIS_mean_2001_2022_kgC_m2 = NPP_MODIS_mean_2001_2022_kgC_m2*10),
                           std.ov = FALSE, 
                           fixed.x = FALSE,
                           group = "aridity_index_factor",
                           group.equal = c("intercepts", "regressions"))


anova(fit_model_SEM, fit_model_SEM_equal) #Multi-group model with some constraints is better than the constrained/equal model
summary(fit_model_SEM, fit.measures = TRUE, standardized = T, rsq = T)


# Example of constraining a path and testing
model_SEM_const0 <- ' 
  POM_mgpg_OC ~ 
      c("clay.poc1", "clay.poc2", "clay.poc3") * clay_total + 
      c("ph.poc1", "ph.poc1", "ph.poc2") * ph_h2o + 
      c("mci.poc1", "mci.poc2", "mci.poc2") * MCI + 
      c("et.poc1", "et.poc2", "et.poc1") * ET + 
      c("map.poc1", "map.poc2", "map.poc1") * prcp_mean +  
      c("npp.poc1", "npp.poc2", "npp.poc2") * NPP_MODIS_mean_2001_2022_kgC_m2 + 
      c("cult.poc1", "cult.poc1", "cult.poc2") * horizon_ap + 
      c("hzn.poc1", "hzn.poc1", "hzn.poc1") * hzn_bot 
    MAOM_mgpg_OC ~ 
      c("clay.maoc1", "clay.maoc1", "clay.maoc1") * clay_total + 
      c("ph.maoc1", "ph.maoc2", "ph.maoc2") * ph_h2o + 
      c("mci.maoc1", "mci.maoc2", "mci.maoc3") * MCI + 
      c("et.maoc1", "et.maoc2", "et.maoc2") * ET +
      c("map.maoc1", "map.maoc2", "map.maoc2") * prcp_mean + 
      c("npp.maoc1", "npp.maoc2", "npp.maoc2") * NPP_MODIS_mean_2001_2022_kgC_m2 + 
      c("cult.maoc1", "cult.maoc1", "cult.maoc1") *  horizon_ap + 
      c("hzn.maoc1", "hzn.maoc1", "hzn.maoc1") * hzn_bot 
    ph_h2o ~ 
      c("et.ph1", "et.ph2", "et.ph2") * ET + 
      c("map.ph1", "map.ph2", "map.ph3") * prcp_mean + 
      c("mci.ph1", "mci.ph2", "mci.ph3") * MCI + 
      c("cult.ph1", "cult.ph2", "cult.ph2") * horizon_ap + 
      c("hzn.ph1", "hzn.ph1", "hzn.ph1") * hzn_bot 
    NPP_MODIS_mean_2001_2022_kgC_m2 ~ 
      c("clay.npp1", "clay.npp2", "clay.npp2") * clay_total +  
      c("ph.npp1", "ph.npp1", "ph.npp1") * ph_h2o + 
      c("et.npp1", "et.npp2", "et.npp3") * ET + 
      c("map.npp1", "map.npp2", "map.npp3") * prcp_mean + 
      c("mci.npp1", "mci.npp2", "mci.npp1") * MCI + 
      c("cult.npp1", "cult.npp2", "cult.npp2") * horizon_ap + 
      c("hzn.npp1", "hzn.npp1", "hzn.npp2") * hzn_bot '

fit_model_SEM_const0 <- sem(model_SEM_const0,
                            data = dat %>%
                              mutate(ET = ET/100,
                                     prcp_mean = prcp_mean/100,
                                     NPP_MODIS_mean_2001_2022_kgC_m2 = NPP_MODIS_mean_2001_2022_kgC_m2*10),
                            std.ov = FALSE, 
                            fixed.x = FALSE,
                            group = "aridity_index_factor")

model_SEM_const1 <- ' 
 POM_mgpg_OC ~ 
      c("clay.poc1", "clay.poc1", "clay.poc1") * clay_total + # Testing this path.
      c("ph.poc1", "ph.poc1", "ph.poc2") * ph_h2o + 
      c("mci.poc1", "mci.poc2", "mci.poc2") * MCI + 
      c("et.poc1", "et.poc2", "et.poc1") * ET + 
      c("map.poc1", "map.poc2", "map.poc1") * prcp_mean +  
      c("npp.poc1", "npp.poc2", "npp.poc2") * NPP_MODIS_mean_2001_2022_kgC_m2 + 
      c("cult.poc1", "cult.poc1", "cult.poc2") * horizon_ap + 
      c("hzn.poc1", "hzn.poc1", "hzn.poc1") * hzn_bot 
    MAOM_mgpg_OC ~ 
      c("clay.maoc1", "clay.maoc1", "clay.maoc1") * clay_total + 
      c("ph.maoc1", "ph.maoc2", "ph.maoc2") * ph_h2o + 
      c("mci.maoc1", "mci.maoc2", "mci.maoc3") * MCI + 
      c("et.maoc1", "et.maoc2", "et.maoc2") * ET +
      c("map.maoc1", "map.maoc2", "map.maoc2") * prcp_mean + 
      c("npp.maoc1", "npp.maoc2", "npp.maoc2") * NPP_MODIS_mean_2001_2022_kgC_m2 + 
      c("cult.maoc1", "cult.maoc1", "cult.maoc1") *  horizon_ap + 
      c("hzn.maoc1", "hzn.maoc1", "hzn.maoc1") * hzn_bot 
    ph_h2o ~ 
      c("et.ph1", "et.ph2", "et.ph2") * ET + 
      c("map.ph1", "map.ph2", "map.ph3") * prcp_mean + 
      c("mci.ph1", "mci.ph2", "mci.ph3") * MCI + 
      c("cult.ph1", "cult.ph2", "cult.ph2") * horizon_ap + 
      c("hzn.ph1", "hzn.ph1", "hzn.ph1") * hzn_bot 
    NPP_MODIS_mean_2001_2022_kgC_m2 ~ 
      c("clay.npp1", "clay.npp2", "clay.npp2") * clay_total +  
      c("ph.npp1", "ph.npp1", "ph.npp1") * ph_h2o + 
      c("et.npp1", "et.npp2", "et.npp3") * ET + 
      c("map.npp1", "map.npp2", "map.npp3") * prcp_mean + 
      c("mci.npp1", "mci.npp2", "mci.npp1") * MCI + 
      c("cult.npp1", "cult.npp2", "cult.npp2") * horizon_ap + 
      c("hzn.npp1", "hzn.npp1", "hzn.npp2") * hzn_bot '

fit_model_SEM_const1 <- sem(model_SEM_const1,
                            data = dat %>%
                              mutate(ET = ET/100,
                                     prcp_mean = prcp_mean/100,
                                     NPP_MODIS_mean_2001_2022_kgC_m2 = NPP_MODIS_mean_2001_2022_kgC_m2*10),
                            std.ov = FALSE, 
                            fixed.x = FALSE,
                            group = "aridity_index_factor")

anova(fit_model_SEM_const0, fit_model_SEM_const1)
#No difference between models when the three clay -> NPP paths are set as equal (i.e.; 1,1,1). 
#If there had been a difference, would have tested for differences between each group (i.e.; 1,2,1; 1,1,3; 1,2,3) 


### Extracting SEM coefficients for figures and tables (Table 1, Table S1) ----
sem_solution_std <- standardizedSolution(fit_model_SEM) 
sem_solution_not_std <- parameterEstimates(fit_model_SEM)

sem_coeffs_std <- sem_solution_std %>%
  filter(pvalue < 0.05 & label != "" & group == "0" & !str_detect(label, ".v.")) %>%
  dplyr::select(lhs, rhs, label, group, est.std) %>%
  mutate(group_new = if_else(group != "0", group,
                             if_else(group == "0" & str_detect(label, "2"), 2, 
                                     if_else(group == "0" & str_detect(label, "3"), 3, 1))),
         Effect = if_else(!str_detect(label, "total") & !str_detect(label, "direct"), "Indirect", 
                          if_else(str_detect(label, "total"), "Total", "Direct")),
         Response = if_else(str_detect(lhs, ".poc"), "POM_mgpg_OC", "MAOM_mgpg_OC"),
         Pred_variable = if_else(Effect == "Direct" & str_detect(label, "ph."), "ph_h2o",
                                 if_else(Effect == "Direct" & str_detect(label, "npp."), "NPP",
                                         if_else(str_detect(label, "cult."), "horizon_ap",
                                                 if_else(str_detect(label, "mci."), "MCI",
                                                         if_else(str_detect(label, "et."), "ET",
                                                                 if_else(str_detect(label, "clay."), "clay_total",
                                                                         if_else(str_detect(label, "map."), "prcp_mean",
                                                                                 if_else(str_detect(label, "hzn."), "hzn_depth",
                                                                                         "ph_h2o")))))))),
         Med_variable = if_else(Effect != "Indirect", "Direct",
                                if_else(Effect == "Indirect" & str_detect(label, ".npp."), "NPP", "pH")),
         group_new = as.factor(group_new)) %>%
  group_by(Effect, Response, Med_variable, Pred_variable, group_new) %>%
  summarize(across(where(is.numeric), list(sum = sum, count = ~ n()))) %>%
  #add blanks for 0s for plotting.
  bind_rows(data.frame(Effect = c(rep("Indirect", 2), rep("Direct", 3)), Med_variable = c(rep("NPP", 2), rep("Direct", 3)), 
                       Pred_variable = c(rep("prcp_mean", 2), rep("NPP", 2), "clay_total"),
                       est.std_sum = 0, group_new = c(rep("2",4), "3"), Response = c(rep(c("POM_mgpg_OC", "MAOM_mgpg_OC"),2), "MAOM_mgpg_OC"))) %>%
  mutate(Pred_variable = factor(Pred_variable, levels = c("horizon_ap", "ET", "prcp_mean", "MCI",
                                                          "clay_total", "NPP", "ph_h2o", "hzn_depth")),
         Pred_variable = fct_recode(Pred_variable, MCI = "MCI", MAP = "prcp_mean",
                                    Cultivation = "horizon_ap", ET = "ET",
                                    pH = "ph_h2o", Texture = "clay_total", Depth = "hzn_depth"),
         Response = factor(Response, levels = c("POM_mgpg_OC", "MAOM_mgpg_OC")),
         Response = fct_recode(Response, POC = "POM_mgpg_OC", MAOC = "MAOM_mgpg_OC")) 

#Figure 2 created outside of R
cat("Fig 2. N = ", summary(fit_model_SEM)$data$nobs %>% sum(), 
    "\n(Omitted", summary(fit_model_SEM)$data$norig %>% sum() - summary(fit_model_SEM)$data$nobs %>% sum(), "due to missing data)")

# 3) Saturation ----
library(quantreg)
library(visreg)

rqfit_1 <- rq(MAOM_mgpg_OC ~ claysilt:aridity_index_factor - 1, 
              tau = 0.95, data = dat)
summary(rqfit_1)
visreg(rqfit_1, "claysilt", by = "aridity_index_factor", overlay=TRUE)


quant_reg_lines_clim <- data.frame(aridity_index_factor = factor(c("Arid", "Sub-humid", "Humid"), levels = c("Arid", "Sub-humid", "Humid")),
                                   intercept = c(0,0,0),
                                   slope = summary(rqfit_1)$coefficients[1:3],
                                   se = summary(rqfit_1)$coefficients[4:6]) %>%
  mutate(ci_lower = slope - se*1.96,
         ci_upper = slope + se*1.96)


MAOM_v_Texture_aridity <-  
  ggplot(dat,
         aes(x = claysilt, y = MAOM_mgpg_OC)) +
  geom_point(aes(color = aridity_index_factor), alpha = 0.5, shape=16, size = 1) +
  geom_point(aes(x = claysilt, y = MAOM_mgpg_OC),
             dat %>% filter(WS_mgpg_OC < 100 & aridity_index_factor == "Arid"),
             alpha = 0.5, color = "#E1BE6A", shape = 16, size = 1) +
  geom_abline(slope = seq(quant_reg_lines_clim$ci_lower[1], 
                          quant_reg_lines_clim$ci_upper[1], 0.001),
              intercept = 0,
              color = "grey80", alpha = .15) +
  geom_abline(slope = seq(quant_reg_lines_clim$ci_lower[2], 
                          quant_reg_lines_clim$ci_upper[2], 0.001),
              intercept = 0,
              color = "grey80", alpha = 0.1) +
  geom_abline(slope = seq(quant_reg_lines_clim$ci_lower[3], 
                          quant_reg_lines_clim$ci_upper[3], 0.001),
              intercept = 0,
              color = "grey80", alpha = 0.1) +
  xlab(expression(paste("Silt and clay content (%)"))) +
  ylab(expression(paste("MAOM C (g C kg"^-1, "soil)")))+
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(limits = c(0, 90), expand = c(0, 0)) +
  geom_abline(aes(slope = slope, intercept = intercept, 
                  color = aridity_index_factor, linetype = aridity_index_factor),
              data = quant_reg_lines_clim, linewidth = .7) +
  theme_classic() +
  scale_fill_manual(name = "Aridity Group",
                    values = c( "#E1BE6A", "#7ec4a6", "#6a9abb")) +
  scale_color_manual(name = "Aridity Group", values = c("#caab5f", "#299d6c", "#064672")) +
  scale_linetype_manual(name = "Aridity Group", values = c(1, 2, 6)) +
  geom_abline(slope = 0.86, intercept = 0, color = "gray60",
              linewidth = .7, linetype = 1) +
  theme(legend.position = c(.2, .85))

cat("Fig 3. N = ", nrow(filter(dat, !is.na(claysilt))), "\n directly measured = ",
    nrow(filter(dat, !is.na(claysilt) & fraction_quant_method == "direct measurement")),
    "MIR predicted = ", nrow(filter(dat, !is.na(claysilt) & fraction_quant_method == "MIR prediction")))

