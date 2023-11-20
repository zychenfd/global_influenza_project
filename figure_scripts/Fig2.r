#==load package==
library(ggplot2)
library(stringr)
library(Biostrings)
library(ape)
library(treeio)
library(adephylo)
library(tidyverse)
library(gtools)
library(dplyr)
library(timeDate)
library(coronavirus)
library(readxl)
library(sf)
library(rgdal)
library(ggtree)
library(ggsci)
library(reshape2)
library(rworldmap)
library(patchwork)
library(ggsci)
library(grid)
library(scales)

#==define color==
colors <- c(pal_npg("nrc", alpha =1)(10)[c(1:7,9:10)],"darkred","#FADDA9","grey80")
colors1 <- c(pal_aaas("default", alpha =0.7)(10))
show_col(colors)
show_col(colors1)
value = c("Japan/Korea" = colors[1],"West Asia" = colors[3],"Western Asia" = colors[3],"Northern America" = colors[6],"North America" = colors[6],
          "Southeast Asia"= colors[4], "South Asia"= colors[5],"South-eastern Asia"= colors[4], "Southern Asia"= colors[5], "Europe"= colors[2], "Oceania"= colors[7],
          "North China" = colors[8], "South China" = colors[11], "Russia"= colors[10],  "Southern America"= colors[12],  "South America"= colors[12], 
          "Africa"= colors[9], "Americas" = colors1[1], "Asia" = colors1[2], "China" = colors1[3])

#==define pathway==
setwd("C:/Users/zyche/Nutstore/1/Evolution_study/Flu_phylogeography/analyses/scripts")

#==read air flow data==
world_data<-getMap(resolution='low')@data
lati_long <- read.csv("C:/Users/zyche/Nutstore/1/Evolution_study/Oxford Visiting/Influenza project/data/lati_long1.csv")

air_data <- read.csv("C:/Users/zyche/Nutstore/1/Evolution_study/Oxford Visiting/Influenza project/data/Air_flow_total.csv") %>%
  filter(region_final_dep != region_final_arr) %>%
  left_join(lati_long, by = c("region_final_dep" = "region_final"))%>%
  left_join(lati_long, by = c("region_final_arr" = "region_final"))

#==read map data==
#1.global map
world_region <- read.csv("C:/Users/zyche/Nutstore/1/Evolution_study/Oxford Visiting/Influenza project/data/world_region.csv")
nine <- st_read("data/Geographic/nine.shp")
nine <- st_transform(nine, CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"))
worldmap0 <- st_read("data/Geographic/World_map.shp")
worldmap <- left_join(worldmap0, world_region[,c(3,4)], by = c("iso_a3" = "alpha.3")) %>% filter(!is.na(region_final))
worldmap1 <- worldmap[,c(16,17)] %>% filter(!region_final %in% c("South China", "China")) %>% group_by(region_final) %>% summarise()
worldmap1 <- st_transform(worldmap1, CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"))
worldmap1$region_final[worldmap1$region_final == "South America"] <- "Southern America"

#2.china map
china_province <- read_xlsx("C:/Users/zyche/Nutstore/1/Evolution_study/Oxford Visiting/Influenza project/data/province.xlsx")
china <- st_read("C:/Users/zyche/Nutstore/1/Evolution_study/Oxford Visiting/Influenza project/data/china_map/gadm41_CHN_1.shp")
china_taiwan <- st_read("C:/Users/zyche/Nutstore/1/Evolution_study/Oxford Visiting/Influenza project/data/china_map/gadm41_TWN_0.shp")
names(china_taiwan)[2] <- "NAME_1"
china <- rbind(china[,c(4,12)], china_taiwan[,c(2,3)])
china$NAME_1 <- str_to_lower(china$NAME_1)
china <- left_join(china, china_province[,1:2], by = c("NAME_1" = "province"))
china$north[china$NAME_1 %in% c("nei mongol", "ningxia hui","xinjiang uygur")] <- "T"
china$north[is.na(china$north)] <- "F"
china1 <- china[,2:3] %>% group_by(north) %>% summarise()
china1 <- st_transform(china1, CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"))
worldmap0 <- st_transform(worldmap0, CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"))
china1$north[china1$north == "T"] <- "North China"
china1$north[china1$north == "F"] <- "South China"

worldmap1$region_final[worldmap1$region_final == "Northern America"] <- "North America"
worldmap1$region_final[worldmap1$region_final == "Southern America"] <- "South America"
worldmap1$region_final[worldmap1$region_final == "South-eastern Asia"] <- "Southeast Asia"
worldmap1$region_final[worldmap1$region_final == "Southern Asia"] <- "South Asia"
worldmap1$region_final[worldmap1$region_final == "Western Asia"] <- "West Asia"
#Supplementary Fig 2
ggplot() +
  geom_sf(data = worldmap0, color = "grey50", fill = "white")+
  geom_sf(data = worldmap1,aes(fill = region_final), color = "black", size = 0.05)+
  geom_sf(data = china1, aes(fill = north), color = "black", size = 0.02)+
  geom_sf(data = nine, color="black",size = 0.01)+
  scale_fill_manual("Geographic regions",values = value)+
  theme_void()+
  theme( legend.text = element_text(size = 7),
         legend.title = element_text(size = 8),
         legend.key.size = unit(0.5,"cm")) -> p_tmp

# svg("output/12_regions.svg",width = 8, height = 3.5)
# p_tmp
# dev.off()
# tiff("output/12_regions.tif",width = 8, height = 3.5,
#      units = "in", res = 326, compression = "lzw")
# p_tmp
# dev.off()
# pdf("output/12_regions.pdf",width = 8, height = 3.5)
# p_tmp
# dev.off()

#==plot==
cut(air_data$Num_per_month, breaks = c(0, 100000, 500000, 1000000, 2000000, 5000000, 10000000), right = T,
    labels = c(0.5,1,2,3,4,5)) ->  air_data$Num_per_month_type
air_data$Num_per_month_type <- as.numeric(as.character(air_data$Num_per_month_type))
ggplot() + 
  geom_sf(data = worldmap0, color = "grey50", fill = "grey90")+
  geom_sf(data = worldmap1,aes(fill = region_final), color = "black", size = 0.1)+
  geom_sf(data = china1, aes(fill = north), color = "black", size = 0.02)+
  geom_sf(data = nine, color="black",size = 0.01)+
  geom_curve(data = air_data[air_data$period == "Period 1" & air_data$Num_per_month >= 100000,],
             aes(x = as.double(long.x), 
                 y = as.double(lati.x), 
                 xend = as.double(long.y), 
                 yend = as.double(lati.y),
                 linewidth = Num_per_month_type),
             color = "#0061A3",
             alpha = 0.5,
             curvature = 0.35)+
  guides(fill = F, linewidth = F)+
  # guides(fill = F)+
  # guides(linewidth = guide_legend(nrow = 1, title.position = "top"))+
  theme_void()+
  geom_point(data = lati_long, aes(x = long, y = lati),shape = 21, size = 3,fill = "lightblue")+
  theme(plot.margin = margin(0,0,0,0, "cm"))+
  scale_x_continuous(limits=c(-170,170))+
  theme(legend.direction = 'horizontal')+
  scale_linewidth_continuous("Average monthly air traffic (million)", range = c(0.3,3), 
                             breaks = seq(1,5,1),
                             labels = c("(0.1, 0.5]","(0.5, 1.0]","(1.0, 2.0]","(2.0, 5.0]","(5.0, 10.0]"), limits = c(1,5))+
  scale_fill_manual(values = value)+
  labs(subtitle = "a. First epoch")-> p1

ggplot() + 
  geom_sf(data = worldmap0, color = "grey50", fill = "grey90")+
  geom_sf(data = worldmap1,aes(fill = region_final), color = "black", size = 0.1)+
  geom_sf(data = china1, aes(fill = north), color = "black", size = 0.02)+
  geom_sf(data = nine, color="black",size = 0.01)+
  geom_curve(data = air_data[air_data$period == "Period 2" & air_data$Num_per_month >= 100000,],
             aes(x = as.double(long.x), 
                 y = as.double(lati.x), 
                 xend = as.double(long.y), 
                 yend = as.double(lati.y),
                 linewidth = Num_per_month_type),
             color = "#0061A3",
             alpha = 0.5,
             curvature = 0.35)+
  # guides(fill = F, linewidth = F)+
  guides(fill = F)+
  # guides(linewidth = guide_legend(nrow = 1, title.position = "top"))+
  theme_void()+
  geom_point(data = lati_long, aes(x = long, y = lati),shape = 21, size = 3,fill = "lightblue")+
  theme(plot.margin = margin(0,0,0,0, "cm"),
        legend.position = "none")+
  scale_x_continuous(limits=c(-170,170))+
  theme(legend.direction = 'horizontal')+
  scale_linewidth_continuous("Average monthly air traffic (million)", range = c(0.3,3), 
                             breaks = seq(1,5,1),
                             labels = c("(0.1, 0.5]","(0.5, 1.0]","(1.0, 2.0]","(2.0, 5.0]","(5.0, 10.0]"), limits = c(1,5))+
  scale_fill_manual(values = value)+
  labs(subtitle = "b. Second epoch") -> p2

ggplot() + 
  geom_sf(data = worldmap0, color = "grey50", fill = "grey90")+
  geom_sf(data = worldmap1,aes(fill = region_final), color = "black", size = 0.1)+
  geom_sf(data = china1, aes(fill = north), color = "black", size = 0.02)+
  geom_sf(data = nine, color="black",size = 0.01)+
  geom_curve(data = air_data[air_data$period == "Period 3" & air_data$Num_per_month >= 100000,],
             aes(x = as.double(long.x), 
                 y = as.double(lati.x), 
                 xend = as.double(long.y), 
                 yend = as.double(lati.y),
                 linewidth = Num_per_month_type),
             color = "#0061A3",
             alpha = 0.5,
             curvature = 0.35)+
  guides(fill = F, linewidth = F)+
  # guides(fill = F)+
  # guides(linewidth = guide_legend(nrow = 1, title.position = "top"))+
  theme_void()+
  geom_point(data = lati_long, aes(x = long, y = lati),shape = 21, size = 3,fill = "lightblue")+
  theme(plot.margin = margin(0,0,0,0, "cm"),
        legend.position = "none")+
  scale_x_continuous(limits=c(-170,170))+
  theme(legend.direction = 'horizontal')+
  scale_linewidth_continuous("Average monthly air traffic (million)", range = c(0.3,3), 
                             breaks = seq(1,5,1),
                             labels = c("(0.1, 0.5]","(0.5, 1.0]","(1.0, 2.0]","(2.0, 5.0]","(5.0, 10.0]"), limits = c(1,5))+
  scale_fill_manual(values = value)+
  labs(subtitle = "c. Third epoch") -> p3

#==Extended Data Figure== 
air_data_s2 <- air_data[air_data$period == "Period 2",]
names(air_data_s2)[6] <- "Num_per_month_2"
air_data_s3 <- air_data[air_data$period == "Period 3",]
names(air_data_s3)[6] <- "Num_per_month_3"
air_data_s <- left_join(air_data[air_data$period == "Period 1",], air_data_s2[,c(2,3,6)])
air_data_s <- left_join(air_data_s, air_data_s3[,c(2,3,6)])
ggplot() + 
  geom_sf(data = worldmap0, color = "grey50", fill = "grey90")+
  geom_sf(data = worldmap1,aes(fill = region_final), color = "black", size = 0.1)+
  geom_sf(data = china1, aes(fill = north), color = "black", size = 0.02)+
  geom_sf(data = nine, color="black",size = 0.01)+
  geom_curve(data = air_data_s,
             aes(x = as.double(long.x), 
                 y = as.double(lati.x), 
                 xend = as.double(long.y), 
                 yend = as.double(lati.y),
                 linewidth = Num_per_month_2/Num_per_month),
             color = "grey",
             alpha = 0.5,
             curvature = 0.35)+
  guides(fill = F)+
  guides(linewidth = guide_legend(nrow = 1, title.position = "top"))+
  theme_void()+
  geom_point(data = lati_long, aes(x = long, y = lati),shape = 21, size = 3,fill = "grey")+
  theme(plot.margin = margin(0,0,0,0, "cm"),
        legend.position = "none")+
  theme(legend.direction = 'horizontal')+
  scale_x_continuous(limits=c(-170,170))+
  scale_fill_manual(values = value)+
  scale_linewidth_continuous("Reduction relative to pre-pandemic level", range = c(0.1,2), 
                             limits = c(0,1.6))+
  labs(subtitle = "a. Second epoch") -> p2_s

ggplot() + 
  geom_sf(data = worldmap0, color = "grey50", fill = "grey90")+
  geom_sf(data = worldmap1,aes(fill = region_final), color = "black", size = 0.1)+
  geom_sf(data = china1, aes(fill = north), color = "black", size = 0.02)+
  geom_sf(data = nine, color="black",size = 0.01)+
  geom_curve(data = air_data_s,
             aes(x = as.double(long.x), 
                 y = as.double(lati.x), 
                 xend = as.double(long.y), 
                 yend = as.double(lati.y),
                 linewidth = Num_per_month_3/Num_per_month),
             color = "grey",
             alpha = 0.5,
             curvature = 0.35)+
  theme_void()+
  geom_point(data = lati_long, aes(x = long, y = lati),shape = 21, size = 3,fill = "grey")+
  theme(plot.margin = margin(0,0,0,0, "cm"),
        legend.position = "none")+
  theme(legend.direction = 'horizontal')+
  scale_x_continuous(limits=c(-170,170))+
  guides(linewidth = guide_legend(nrow = 1, title.position = "top"))+
  guides(fill = F, linewidth = F)+
  scale_fill_manual(values = value)+
  scale_linewidth_continuous("Reduction relative to pre-pandemic level", range = c(0.1,2), 
                             limits = c(0,1.6))+
  labs(subtitle = "b. Third epoch") -> p3_s

pdf("output/Air_reduction.pdf",width = 7.5, height = 3)
((p2_s + p3_s) + plot_layout(nrow = 1,ncol = 2,heights = c(1,1),guides = "collect")&theme(legend.position = "bottom"))
dev.off()

#==Air traffic comparison between each pair of regions across epochs
air_data1 <- air_data 
for (i in 1:nrow(air_data1)) {
  air_data1$route[i] <- paste0(sort(c(air_data1$region_final_dep[i],air_data1$region_final_arr[i])),collapse = "_")
}

air_data1 <- air_data1 %>%
  group_by(period, route) %>%
  summarise(dual_no  = sum(Num_per_month))
air_data2 <- reshape2::dcast(air_data1, route~period)
air_data2$continent_a <- sapply(str_split(air_data2$route,"_"), function(x) x[1])
air_data2$continent_b <- sapply(str_split(air_data2$route,"_"), function(x) x[2])
air_data2$continent_a[air_data2$continent_a %in% c("Russia")] <- "Europe"
air_data2$continent_a[air_data2$continent_a %in% c("Northern America", "South America")] <- "Americas"
air_data2$continent_a[air_data2$continent_a %in% c("South-eastern Asia","Southern Asia","Japan/Korea","Western Asia","South China", "North China")] <- "Asia"
air_data2$continent_b[air_data2$continent_b %in% c("Russia")] <- "Europe"
air_data2$continent_b[air_data2$continent_b %in% c("Northern America", "South America")] <- "Americas"
air_data2$continent_b[air_data2$continent_b %in% c("South-eastern Asia","Southern Asia","Japan/Korea","Western Asia","South China", "North China")] <- "Asia"

for (i in 1:nrow(air_data2)) {
  air_data2$route1[i] <- paste0(sort(c(air_data2$continent_a[i],air_data2$continent_b[i])),collapse = " - ")
}

min(log(air_data2$`Period 1`))
max(log(air_data2$`Period 1`))
min(log(air_data2$`Period 2`))
max(log(air_data2$`Period 2`))
min(log(air_data2$`Period 3`))
max(log(air_data2$`Period 3`))

air_data2$Period_1 <- log(air_data2$`Period 1`)
air_data2$Period_2 <- log(air_data2$`Period 2`)
air_data2$Period_3 <- log(air_data2$`Period 3`)
air_data2$prop_2_1 <- air_data2$`Period 2`/air_data2$`Period 1`
air_data2$prop_3_1 <- air_data2$`Period 3`/air_data2$`Period 1`
range(air_data2$prop_2_1)
range(air_data2$prop_3_1)
air_data2$route_type[air_data2$continent_a != air_data2$continent_b] <- "Inter-continental" 
air_data2$route_type[air_data2$continent_a == air_data2$continent_b] <- "Intra-continental"
air_data2$whe_rec2_1 <- "No"
air_data2$whe_rec2_1[air_data2$prop_2_1 > 1] <- "Yes"
air_data2$whe_rec3_1 <- "No"
air_data2$whe_rec3_1[air_data2$prop_3_1 > 1] <- "Yes"

ggplot(data = air_data2) +
  annotate("segment", x = 6, xend = 17.5, y = 6, yend = 17.5, linetype = 2, size = 0.2)+
  annotate("segment", x = 6, xend = 16.5, y = 7, yend = 17.5, linetype = 2, size = 0.2)+
  annotate("segment", x = 6, xend = 15.5, y = 8, yend = 17.5, linetype = 2, size = 0.2)+
  annotate("segment", x = 6, xend = 14.5, y = 9, yend = 17.5, linetype = 2, size = 0.2)+
  annotate("segment", x = log(14868145.4), xend = log(14868145.4), y = log(18682667.0), yend = 14.5,  size = 0.2)+
  annotate("segment", x = log(3833.6), xend = log(3833.6), y = log(115832.9), yend = 12.5,  size = 0.2)+
  geom_point(aes(x = log(`Period 2`), y = log(`Period 1`), fill = route_type, size = prop_2_1*100), shape = 21)+
  # geom_point(data = air_data2[air_data2$whe_rec2_1 == "Yes",], aes(x = log(`Period 2`), y = log(`Period 1`), fill = route_type, size = prop_2_1*100), shape = 23)+
  scale_x_continuous("Epoch 2", limits = c(6,17.5), expand = c(0,0))+
  scale_y_continuous("Epoch 1", limits = c(6,17.5), expand = c(0,0))+
  scale_size_continuous("Proportional changes\nrelative to epoch 1 (%)",
                        limits = c(min(air_data2$prop_2_1 * 100), max(air_data2$prop_3_1 * 100)))+
  scale_fill_manual("Air traffic types", values = colors1[c(4,5)])+
  # scale_shape_binned(values = c(21,23))+
  theme_bw()+
  guides(fill = guide_legend(ncol = 1),
         size = guide_legend(nrow = 1))+
  annotate("text", x= 15.7, y = 14.5-0.1, label = "North China -\nSouth China", size = 2, vjust = 1)+
  annotate("text", x= log(3833.6), y = 12.5+0.1, label = "North China -\nRussia", size = 2, vjust = 0)+
  theme(panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent"),
        legend.spacing.y = unit(0,"char"),
        plot.background =  element_rect(fill = "transparent"))+
  labs(subtitle = "Log-trans air traffic volume",tag = "e")-> c1

ggplot(data = air_data2) +
  annotate("segment", x = 6, xend = 17.5, y = 6, yend = 17.5, linetype = 2, size = 0.2)+
  annotate("segment", x = 6, xend = 16.5, y = 7, yend = 17.5, linetype = 2, size = 0.2)+
  annotate("segment", x = 6, xend = 15.5, y = 8, yend = 17.5, linetype = 2, size = 0.2)+
  annotate("segment", x = 6, xend = 14.5, y = 9, yend = 17.5, linetype = 2, size = 0.2)+
  annotate("segment", x = 12, xend = 13, y = log(59106.5), yend = log(59106.5),  size = 0.2)+
  annotate("segment", x = log(25818.5), xend = log(25818.5), y = log(324891.8), yend = 14.5,  size = 0.2)+
  annotate("segment", x = log(65395.4), xend = log(65395.4), y = log(63695.4), yend = 8.6,  size = 0.2)+
  geom_point(aes(x = log(`Period 3`), y = log(`Period 1`), fill = route_type, size = prop_3_1*100), shape =21)+
  geom_point(data = air_data2[air_data2$whe_rec3_1 == "Yes",],aes(x = log(`Period 3`), y = log(`Period 1`), size = prop_3_1*100), shape =21, fill = "transparent", color = "red")+
  # geom_point(data = air_data2[air_data2$whe_rec3_1 == "Yes",], aes(x = log(`Period 3`), y = log(`Period 1`), fill = route_type, size = prop_2_1*100), shape = 24)+
  scale_x_continuous("Epoch 3", limits = c(6,17.5), expand = c(0,0))+
  scale_y_continuous("Epoch 1", limits = c(6,17.5), expand = c(0,0))+
  scale_size_continuous("Proportional changes\nrelative to epoch 1 (%)",
                        limits = c(min(air_data2$prop_2_1 * 100), max(air_data2$prop_3_1 * 100)))+
  scale_fill_manual("Air traffic types", values = colors1[c(4,5)])+
  theme_bw()+
  guides(fill = F,
         size = F)+
  annotate("text", x= 13.1, y = log(59106.5), label = "Africa -\nRussia", size = 2, vjust = 0.5, hjust = 0)+
  annotate("text", x= log(65395.4), y = 8.5, label = "South America -\nWestern Asia", size = 2, vjust = 1)+
  annotate("text", x= log(25818.5), y = 14.5+0.1, label = "North China -\nNorthern America", size = 2, vjust = 0)+
  theme(panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent"),
        legend.spacing.y = unit(0,"char"),
        plot.background =  element_rect(fill = "transparent"))+
  labs(subtitle = "Log-trans air traffic volume",tag = "f")-> c2


#==air traffic from and to each region==
list.files("C:/Users/zyche/Nutstore/1/OAG_DATA/AIRPOT2AIRPOT", full.names = T) -> tmp_files
tmp_files <- tmp_files[c(1:(length(tmp_files)-1))]

#CHINA airport
china_province <- read_xlsx("C:/Users/zyche/Nutstore/1/Evolution_study/Oxford Visiting/Influenza project/data/province.xlsx")
china_air <- read_xlsx("C:/Users/zyche/Nutstore/1/Evolution_study/Oxford Visiting/Influenza project/data/China_airport.xlsx") %>% 
  filter(!is.na(Province)) %>%
  dplyr::select(c("Province","City","IATA")) %>%
  left_join(china_province, by = c("Province" = "province"))

world_region$region_final_dep <- world_region$region_final 
world_region$region_final_arr <- world_region$region_final 

#==4.1 air flow between regions==
library(data.table)
n_volume <- c()
for (i in 1:length(tmp_files)) {
  # i = 3
  (tmp_files)[i] -> tmp_path
  fread(tmp_path) -> tmp_data
  data.frame(tmp_data) -> tmp_data
  
  tmp_data <- tmp_data %>%
    filter(Dep.Country.Code %in% world_region$alpha.2) %>%
    filter(Arr.Country.Code %in% world_region$alpha.2) %>%
    left_join(world_region[,c(2,5)], by = c("Dep.Country.Code" = "alpha.2")) %>%
    left_join(world_region[,c(2,6)], by = c("Arr.Country.Code" = "alpha.2"))
  
  tmp_data$region_final_dep[tmp_data$Dep.Airport.Code %in% china_air$IATA[china_air$north == "T"] &
                              tmp_data$Dep.Country.Code == "CN"] <- "North China"
  tmp_data$region_final_dep[tmp_data$Dep.Airport.Code %in% china_air$IATA[china_air$north == "F"] &
                              tmp_data$Dep.Country.Code == "CN"] <- "South China"
  tmp_data$region_final_arr[tmp_data$Arr.Airport.Code %in% china_air$IATA[china_air$north == "T"] &
                              tmp_data$Arr.Country.Code == "CN"] <- "North China"
  tmp_data$region_final_arr[tmp_data$Arr.Airport.Code %in% china_air$IATA[china_air$north == "F"] &
                              tmp_data$Arr.Country.Code == "CN"] <- "South China"
  
  tmp_data <- tmp_data %>%
    filter(!(region_final_dep == "South China" & region_final_arr == "North China")) %>%
    filter(!(region_final_dep == "North China" & region_final_arr == "South China")) # remove air flow within China (between North China and South China)
  
  tmp_data1 <- tmp_data %>%
    filter(region_final_dep != region_final_arr) %>%
    group_by(Time.Series,region_final_dep,region_final_arr) %>%
    summarise(Volume = sum(Total.Est..Pax))
  
  bind_rows(n_volume, tmp_data1) -> n_volume
}

#plot
n_volume1 <- n_volume %>% group_by(Time.Series, region_final_dep) %>% summarise(Volume_dep = sum(Volume))
n_volume2 <- n_volume %>% group_by(Time.Series, region_final_arr) %>% summarise(Volume_arr = sum(Volume))
n_volume3 <- left_join(n_volume1, n_volume2, by = c("Time.Series" = "Time.Series","region_final_dep" = "region_final_arr"))
n_volume3$Volume_total <- n_volume3$Volume_dep + n_volume3$Volume_arr
# write.csv(n_volume3, "../data/air_flow.csv", row.names = F)
n_volume3$date <- as.Date(paste0(n_volume3$Time.Series,"15"), format = "%Y%m%d")
n_volume3$region_final_dep[n_volume3$region_final_dep == "South America"] <- "Southern America"
ggplot(data = n_volume3) +
  annotate("rect", xmin = as.Date("2020-02-01"),xmax = as.Date("2021-07-31"),
           ymin = -1000000,ymax = 32000000,alpha = 0.1,fill = "red")+
  geom_line(aes(x = date, y = Volume_total, color = region_final_dep))+
  theme_bw()+
  scale_x_date(date_breaks = "1 year", date_labels = "%Y",
               limits = c(as.Date("2018-01-01"),as.Date("2023-08-01")))+
  labs(x = "Date", y = "Air traffic volumes (m)", subtitle = "Air traffic from and to that region", tag = "d")+
  scale_y_continuous(labels = seq(0,30,10), limits = c(-1000000,32000000), expand = c(0,0))+
  scale_color_manual("Geographic regions",values = value)+
  guides(color = F)+
  theme(legend.position = "none",
        # axis.text.x = element_text(angle = 45, hjust = 1 ,vjust = 1),
        panel.grid.minor = element_blank())-> c0

#==for stat analyses
library(data.table)
n_volume <- c()
for (i in 1:length(tmp_files)) {
  # i = 3
  (tmp_files)[i] -> tmp_path
  fread(tmp_path) -> tmp_data
  data.frame(tmp_data) -> tmp_data
  
  tmp_data <- tmp_data %>%
    filter(Dep.Country.Code %in% world_region$alpha.2) %>%
    filter(Arr.Country.Code %in% world_region$alpha.2) %>%
    left_join(world_region[,c(2,5)], by = c("Dep.Country.Code" = "alpha.2")) %>%
    left_join(world_region[,c(2,6)], by = c("Arr.Country.Code" = "alpha.2"))
  
  # unique(tmp_data$region_final_dep)
  # unique(tmp_data$region_final_arr)
  tmp_data$region_final_dep[str_detect(tmp_data$region_final_dep, "China|Asia|Japan")] <- "Asia"
  tmp_data$region_final_arr[str_detect(tmp_data$region_final_arr, "China|Asia|Japan")] <- "Asia"
  
  tmp_data1 <- tmp_data %>%
    filter(region_final_dep != region_final_arr) %>%
    group_by(Time.Series,region_final_dep,region_final_arr) %>%
    summarise(Volume = sum(Total.Est..Pax))
  
  bind_rows(n_volume, tmp_data1) -> n_volume
}
n_volume1 <- n_volume %>% group_by(Time.Series, region_final_dep) %>% summarise(Volume_dep = sum(Volume))
n_volume2 <- n_volume %>% group_by(Time.Series, region_final_arr) %>% summarise(Volume_arr = sum(Volume))
n_volume3 <- left_join(n_volume1, n_volume2, by = c("Time.Series" = "Time.Series","region_final_dep" = "region_final_arr"))
n_volume3$Volume_total <- n_volume3$Volume_dep + n_volume3$Volume_arr

#==Predictors output from GLM phylogeography==
text <- c("Sample size (D)", "Sample size (O)", "Infect dis vulner index (D)", "Infect dis vulner index (O)", 
          "Stringency index (D)", "Stringency index (O)",  "Internal air traffic (D)", "Internal air traffic (O)",  "Air traffic", 
          "Relative humidity (D)", "Relative humidity (O)", "Precipitation (D)",  "Precipitation (O)", 
          "Temperature (D)", "Temperature (O)", "Latitude (D)", "Latitude (O)", "Distance",                
          "Prop age under 15 (D)", "Prop age under 15 (O)", "Pop density (D)",  "Pop density (O)", "Pop size (D)", "Pop size (O)")

#==
h1n1_even <- data.frame(rbind(t(read.delim("../phylogeography/4.1post-analyses/glm_log/h1n1_even_glm_log.tsv"))))
h1n1_even$type <- row.names(h1n1_even)
h1n1_even <- h1n1_even[,c(1,8,12)]
colnames(h1n1_even) <- h1n1_even[1,]
h1n1_even$type <- "H1N1pdm09"
h3n2_even <- data.frame(rbind(t(read.delim("../phylogeography/4.1post-analyses/glm_log/h3n2_even_glm_log.tsv"))))
h3n2_even$type <- row.names(h3n2_even)
h3n2_even <- h3n2_even[,c(1,8,12)]
colnames(h3n2_even) <- h3n2_even[1,]
h3n2_even$type <- "H3N2"
bv_even <- data.frame(rbind(t(read.delim("../phylogeography/4.1post-analyses/glm_log/bv_even_glm_log.tsv"))))
bv_even$type <- row.names(bv_even)
bv_even <- bv_even[,c(1,8,12)]
colnames(bv_even) <- bv_even[1,]
bv_even$type <- "B/Victoria"
by_even <- data.frame(rbind(t(read.delim("../phylogeography/4.1post-analyses/glm_log/by_even_glm_log.tsv"))))
by_even$type <- row.names(by_even)
by_even <- by_even[,c(1,8,12)]
colnames(by_even) <- by_even[1,]
by_even$type <- "B/Yamagata"

even <- rbind(h1n1_even, h3n2_even, bv_even, by_even)
rownames(even) <- NULL
even <- even %>% 
  filter(mean != "mean") %>%
  filter(str_detect(Summary.Statistic,"Times")) %>%
  mutate(mean = as.numeric(mean)) %>%
  mutate(`95% HPD interval` = str_remove_all(`95% HPD interval`,"\\[|\\]")) %>%
  mutate(HPD_low = as.numeric(str_trim(sapply(str_split(`95% HPD interval`,","), function(x) x[1])))) %>%
  mutate(HPD_upp = as.numeric(str_trim(sapply(str_split(`95% HPD interval`,","), function(x) x[2])))) %>%
  mutate(predictor = rep(1:24,4))

factor(even$predictor, levels = 1:24) -> even$predictor
factor(even$type, levels = unique(even$type)) -> even$type
ggplot(data = even) +
  annotate("rect", xmin = 15.5,xmax = 16.5,
           ymin = -0.5,ymax = 1.5,fill = "lightblue")+
  geom_hline(yintercept = 0, linetype = 2, size = 0.2)+
  geom_errorbar(aes(x = predictor, ymin = HPD_low, ymax = HPD_upp, group = type), width = 0.5, position = position_dodge(width = 0.7), size = 0.3)+
  geom_point(aes(x = predictor, y = mean, fill = type), shape = 21, size = 2,alpha = 1, position = position_dodge(width = 0.7))+
  # coord_flip(clip = "off")+
  scale_x_discrete(labels = rev(text))+
  scale_y_continuous(limits = c(-0.5,1.5),expand = c(0,0))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.title = element_blank(),
        panel.background = element_rect(fill = "transparent",color = "transparent"),
        panel.border = element_rect(fill = "transparent"),
        plot.background =  element_rect(fill = "transparent",color = "transparent"),
        legend.background = element_rect(fill = "transparent",color = "transparent"),
        legend.position = c(0.1, 0.75))+
  labs(x = "", y = "Coefficient * Inclusion probability",tag = "g ")+
  scale_fill_manual("",values = colors1[c(1:3,9)])-> p4

#==figure output==
pdf("output/Fig2.pdf",width = 10, height = 8)
((p1 + p2 + p3) + plot_layout(nrow = 1,ncol = 3,heights = c(1,1,1))&theme(legend.position = "bottom")) -> part1
((c0 + c1 + c2) + plot_layout(nrow = 1,ncol = 3,widths = c(1, 0.7,0.7), guides = "collect")&theme(legend.position = "right")) -> part2
p4 -> part3
viewport(x = 0.005, y = 0.71, width = 1, height = 0.29, just = c("left", "bottom")) -> vp1
viewport(x = 0, y = 0.39, width = 1, height = 0.34, just = c("left", "bottom")) -> vp2
viewport(x = 0, y = -0.03, width = 1, height = 0.45, just = c("left", "bottom")) -> vp3
print(part2,vp = vp2)
print(part1,vp = vp1)
print(part3,vp = vp3)
dev.off()

tiff("output/Fig2.tif",width = 10, height = 8, units = "in", res = 400, compression = "lzw")
((p1 + p2 + p3) + plot_layout(nrow = 1,ncol = 3,heights = c(1,1,1))&theme(legend.position = "bottom")) -> part1
((c0 + c1 + c2) + plot_layout(nrow = 1,ncol = 3,widths = c(1, 0.7,0.7), guides = "collect")&theme(legend.position = "right")) -> part2
p4 -> part3
viewport(x = 0.005, y = 0.71, width = 1, height = 0.29, just = c("left", "bottom")) -> vp1
viewport(x = 0, y = 0.39, width = 1, height = 0.34, just = c("left", "bottom")) -> vp2
viewport(x = 0, y = -0.03, width = 1, height = 0.45, just = c("left", "bottom")) -> vp3
print(part2,vp = vp2)
print(part1,vp = vp1)
print(part3,vp = vp3)
dev.off()

svg("output/Fig2.svg",width = 10, height = 8)
((p1 + p2 + p3) + plot_layout(nrow = 1,ncol = 3,heights = c(1,1,1))&theme(legend.position = "bottom")) -> part1
((c0 + c1 + c2) + plot_layout(nrow = 1,ncol = 3,widths = c(1, 0.7,0.7), guides = "collect")&theme(legend.position = "right")) -> part2
p4 -> part3
viewport(x = 0.005, y = 0.71, width = 1, height = 0.29, just = c("left", "bottom")) -> vp1
viewport(x = 0, y = 0.39, width = 1, height = 0.34, just = c("left", "bottom")) -> vp2
viewport(x = 0, y = -0.03, width = 1, height = 0.45, just = c("left", "bottom")) -> vp3
print(part2,vp = vp2)
print(part1,vp = vp1)
print(part3,vp = vp3)
dev.off()