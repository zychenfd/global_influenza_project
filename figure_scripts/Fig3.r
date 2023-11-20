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
value = c("9.Japan/Korea" = colors[1],"7.WesternAsia" = colors[3],"4.NorthernAmerica" = colors[6],
          "93.SoutheasternAsia"= colors[4], "8.SouthernAsia"= colors[5], "2.Europe"= colors[2], "6.Oceania"= colors[7],
          "91.NorthChina" = colors[8], "92.SouthChina" = colors[11], "3.Russia"= colors[10],  "5.SouthernAmerica"= colors[12], 
          "1.Africa"= colors[9], "Americas" = colors1[1], "Asia" = colors1[2], "China" = colors1[3])
label_region <- c("Africa","Europe","Russia","North America","South America",
                  "Oceania", "West Asia", "South Asia", "Japan/Korea", "North China",
                  "South China", "Southeast Asia")

#==define pathway==
setwd("C:/Users/zyche/Nutstore/1/Evolution_study/Flu_phylogeography/analyses/scripts")

#==read data and plot==
trunk_h1n1_even <- read.delim("../phylogeography/4.1post-analyses/trunk/out_h1n1_even_1.skylines") %>% mutate(type = "H1N1")
trunk_h3n2_even <- read.delim("../phylogeography/4.1post-analyses/trunk/out_h3n2_even_1.skylines") %>% mutate(type = "H3N2")
trunk_bv_even <- read.delim("../phylogeography/4.1post-analyses/trunk/out_bv_even_1.skylines") %>% mutate(type = "BV")
trunk_by_even <- read.delim("../phylogeography/4.1post-analyses/trunk/out_by_even_1.skylines") %>% mutate(type = "BY")

trunk_even <- rbind(trunk_h1n1_even, trunk_h3n2_even, trunk_bv_even, trunk_by_even)
trunk_even$statistic <- str_replace_all(trunk_even$statistic,"pro_","")
trunk_even$date <- decimal2Date(trunk_even$time)

trunk_h1n1_tempo <- read.delim("../phylogeography/4.1post-analyses/trunk/out_h1n1_tempo_1.skylines") %>% mutate(type = "H1N1")
trunk_h3n2_tempo <- read.delim("../phylogeography/4.1post-analyses/trunk/out_h3n2_tempo_1.skylines") %>% mutate(type = "H3N2")
trunk_bv_tempo <- read.delim("../phylogeography/4.1post-analyses/trunk/out_bv_tempo_1.skylines") %>% mutate(type = "BV")
trunk_by_tempo <- read.delim("../phylogeography/4.1post-analyses/trunk/out_by_tempo_1.skylines") %>% mutate(type = "BY")

trunk_tempo <- rbind(trunk_h1n1_tempo, trunk_h3n2_tempo, trunk_bv_tempo, trunk_by_tempo)
trunk_tempo$statistic <- str_replace_all(trunk_tempo$statistic,"pro_","")
trunk_tempo$date <- decimal2Date(trunk_tempo$time)

trunk_even <- left_join(trunk_even, trunk_tempo[,c(1,2,4,6)], by = c("statistic" = "statistic", "time" = "time", "type" = "type"))
trunk_even$mean <- (trunk_even$mean.x + trunk_even$mean.y)/2

trunk_even1 <- trunk_even %>%
  mutate(mean1 = ifelse(mean > 0, mean, 0.000001)) %>%
  group_by(type, date) %>%
  summarise(Shannon_index = sum(- mean1 * log(mean1)))

# trunk_even$statistic[trunk_even$statistic == "NorthChina"] <- "North China"
# trunk_even$statistic[trunk_even$statistic == "NorthernAmerica"] <- "Northern America"
# trunk_even$statistic[trunk_even$statistic == "SouthChina"] <- "South China"
# trunk_even$statistic[trunk_even$statistic == "SoutheasternAsia"] <- "South-eastern Asia"
# trunk_even$statistic[trunk_even$statistic == "SouthernAmerica"] <- "Southern America"
# trunk_even$statistic[trunk_even$statistic == "SouthernAsia"] <- "Southern Asia"
# trunk_even$statistic[trunk_even$statistic == "WesternAsia"] <- "Western Asia"
factor(trunk_even$statistic, levels = unique(trunk_even$statistic)[c(1,2,7,5,10,6,12,11,3,4,8,9)]) -> trunk_even$statistic
# trunk_even <- trunk_even[,c(1,2,4,6)]
trunk_even <- trunk_even[,c(1,2,9,6)]

trunk_even1_h1n1 <- dcast(trunk_even[trunk_even$type == "H1N1",], time ~ statistic, value.var = "mean")
trunk_even1_h1n1$date <- decimal2Date(trunk_even1_h1n1$time)
names(trunk_even1_h1n1)[10] <- "JapanKorea"
ggplot(trunk_even1_h1n1)+
  geom_ribbon(aes(x = date, ymin = 0, ymax = 1, fill = "1.Africa"))+
  geom_ribbon(aes(x = date, ymin = 1-Africa-Europe, ymax = 1-Africa, fill = "2.Europe"))+
  geom_ribbon(aes(x = date, ymin = 1-Africa-Europe-Russia, ymax = 1-Africa-Europe, fill = "3.Russia"))+
  geom_ribbon(aes(x = date, ymin = 1-Africa-Europe-Russia-NorthernAmerica, 
                  ymax = 1-Africa-Europe-Russia, fill = "4.NorthernAmerica"))+
  geom_ribbon(aes(x = date, ymin = 1-Africa-Europe-Russia-NorthernAmerica-SouthernAmerica, 
                  ymax = 1-Africa-Europe-Russia-NorthernAmerica, fill = "5.SouthernAmerica"))+
  geom_ribbon(aes(x = date, ymin = 1-Africa-Europe-Russia-NorthernAmerica-SouthernAmerica-Oceania, 
                  ymax = 1-Africa-Europe-Russia-NorthernAmerica-SouthernAmerica, fill = "6.Oceania"))+
  geom_ribbon(aes(x = date, ymin = 1-Africa-Europe-Russia-NorthernAmerica-SouthernAmerica-Oceania-WesternAsia, 
                  ymax = 1-Africa-Europe-Russia-NorthernAmerica-SouthernAmerica-Oceania, fill = "7.WesternAsia"))+
  geom_ribbon(aes(x = date, ymin = 1-Africa-Europe-Russia-NorthernAmerica-SouthernAmerica-Oceania-WesternAsia-SouthernAsia, 
                  ymax = 1-Africa-Europe-Russia-NorthernAmerica-SouthernAmerica-Oceania-WesternAsia, fill = "8.SouthernAsia"))+
  geom_ribbon(aes(x = date, ymin = NorthChina+SouthChina+SoutheasternAsia, 
                  ymax = JapanKorea+NorthChina+SouthChina+SoutheasternAsia, fill = "9.Japan/Korea"))+
  geom_ribbon(aes(x = date, ymin = SouthChina+SoutheasternAsia, 
                  ymax = NorthChina+SouthChina+SoutheasternAsia, fill = "91.NorthChina"))+
  geom_ribbon(aes(x = date, ymin = SoutheasternAsia, 
                  ymax = SouthChina+SoutheasternAsia, fill = "92.SouthChina"))+
  geom_ribbon(aes(x = date, ymin = 0, 
                  ymax = SoutheasternAsia, fill = "93.SoutheasternAsia"))+
  geom_vline(xintercept = c(as.Date("2020-02-01"),as.Date("2021-08-01")), linetype = 2, linewidth = 0.2)+
  geom_line(data = trunk_even1[trunk_even1$type == "H1N1",],
            aes(x = date, y = Shannon_index/2.5), linetype = 4, color = "grey20")+
  scale_fill_manual("",values = value, label = label_region)+
  scale_x_date("Date", date_breaks = "6 month",date_labels = "%b %Y",
               limits = c(as.Date("2018-02-01"),as.Date("2022-08-01")), expand = c(0.01,0))+
  annotate("text", x = c(as.Date("2019-02-15"),as.Date("2020-11-01"),as.Date("2022-02-01")), y = 1.1, label = c("First epoch","Second epoch","Third epoch"))+
  annotate("segment", x = c(as.Date("2020-02-01"),as.Date("2021-08-01")),xend = c(as.Date("2020-02-01"),as.Date("2021-08-01")),
           yend = 1.02, y = 1.17, arrow = arrow(type = "open", length = unit(0.05, "inches")))+
  coord_cartesian(ylim = c(0,1), clip = "off")+
  scale_y_continuous(expand = c(0.02,0),sec.axis = sec_axis(trans=~.* 2.5,name="Shannon index"))+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        panel.border = element_rect(linewidth = 0.3),
        panel.grid = element_blank())+
  # guides(fill = guide_legend(ncol = 6, nrow = 2))+
  labs(subtitle = "H1N1pdm09\n", x = "",y = "Trunk probablity")-> p1

trunk_even1_h3n2 <- dcast(trunk_even[trunk_even$type == "H3N2",], time ~ statistic, value.var = "mean")
trunk_even1_h3n2$date <- decimal2Date(trunk_even1_h3n2$time)
names(trunk_even1_h3n2)[10] <- "JapanKorea"
ggplot(trunk_even1_h3n2)+
  geom_ribbon(aes(x = date, ymin = 0, ymax = 1, fill = "1.Africa"))+
  geom_ribbon(aes(x = date, ymin = 1-Africa-Europe, ymax = 1-Africa, fill = "2.Europe"))+
  geom_ribbon(aes(x = date, ymin = 1-Africa-Europe-Russia, ymax = 1-Africa-Europe, fill = "3.Russia"))+
  geom_ribbon(aes(x = date, ymin = 1-Africa-Europe-Russia-NorthernAmerica, 
                  ymax = 1-Africa-Europe-Russia, fill = "4.NorthernAmerica"))+
  geom_ribbon(aes(x = date, ymin = 1-Africa-Europe-Russia-NorthernAmerica-SouthernAmerica, 
                  ymax = 1-Africa-Europe-Russia-NorthernAmerica, fill = "5.SouthernAmerica"))+
  geom_ribbon(aes(x = date, ymin = 1-Africa-Europe-Russia-NorthernAmerica-SouthernAmerica-Oceania, 
                  ymax = 1-Africa-Europe-Russia-NorthernAmerica-SouthernAmerica, fill = "6.Oceania"))+
  geom_ribbon(aes(x = date, ymin = 1-Africa-Europe-Russia-NorthernAmerica-SouthernAmerica-Oceania-WesternAsia, 
                  ymax = 1-Africa-Europe-Russia-NorthernAmerica-SouthernAmerica-Oceania, fill = "7.WesternAsia"))+
  geom_ribbon(aes(x = date, ymin = 1-Africa-Europe-Russia-NorthernAmerica-SouthernAmerica-Oceania-WesternAsia-SouthernAsia, 
                  ymax = 1-Africa-Europe-Russia-NorthernAmerica-SouthernAmerica-Oceania-WesternAsia, fill = "8.SouthernAsia"))+
  geom_ribbon(aes(x = date, ymin = NorthChina+SouthChina+SoutheasternAsia, 
                  ymax = JapanKorea+NorthChina+SouthChina+SoutheasternAsia, fill = "9.Japan/Korea"))+
  geom_ribbon(aes(x = date, ymin = SouthChina+SoutheasternAsia, 
                  ymax = NorthChina+SouthChina+SoutheasternAsia, fill = "91.NorthChina"))+
  geom_ribbon(aes(x = date, ymin = SoutheasternAsia, 
                  ymax = SouthChina+SoutheasternAsia, fill = "92.SouthChina"))+
  geom_ribbon(aes(x = date, ymin = 0, 
                  ymax = SoutheasternAsia, fill = "93.SoutheasternAsia"))+
  geom_vline(xintercept = c(as.Date("2020-02-01"),as.Date("2021-08-01")), linetype = 2, linewidth = 0.2)+
  geom_line(data = trunk_even1[trunk_even1$type == "H3N2",],
            aes(x = date, y = Shannon_index/2.5), linetype = 4, color = "grey20")+
  scale_fill_manual("",values = value, label = label_region)+
  scale_x_date("Date", date_breaks = "6 month",date_labels = "%b %Y",
               limits = c(as.Date("2018-02-01"),as.Date("2022-08-01")), expand = c(0.01,0))+
  scale_y_continuous(expand = c(0.02,0),sec.axis = sec_axis(trans=~.* 2.5,name="Shannon index"))+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        panel.grid = element_blank())+
  # guides(fill = guide_legend(ncol = 6, nrow = 2))+
  labs(subtitle = "H3N2",x = "",y = "Trunk probablity")-> p2

trunk_even1_bv <- dcast(trunk_even[trunk_even$type == "BV",], time ~ statistic, value.var = "mean")
trunk_even1_bv$date <- decimal2Date(trunk_even1_bv$time)
names(trunk_even1_bv)[10] <- "JapanKorea"
ggplot(trunk_even1_bv)+
  geom_ribbon(aes(x = date, ymin = 0, ymax = 1, fill = "1.Africa"))+
  geom_ribbon(aes(x = date, ymin = 1-Africa-Europe, ymax = 1-Africa, fill = "2.Europe"))+
  geom_ribbon(aes(x = date, ymin = 1-Africa-Europe-Russia, ymax = 1-Africa-Europe, fill = "3.Russia"))+
  geom_ribbon(aes(x = date, ymin = 1-Africa-Europe-Russia-NorthernAmerica, 
                  ymax = 1-Africa-Europe-Russia, fill = "4.NorthernAmerica"))+
  geom_ribbon(aes(x = date, ymin = 1-Africa-Europe-Russia-NorthernAmerica-SouthernAmerica, 
                  ymax = 1-Africa-Europe-Russia-NorthernAmerica, fill = "5.SouthernAmerica"))+
  geom_ribbon(aes(x = date, ymin = 1-Africa-Europe-Russia-NorthernAmerica-SouthernAmerica-Oceania, 
                  ymax = 1-Africa-Europe-Russia-NorthernAmerica-SouthernAmerica, fill = "6.Oceania"))+
  geom_ribbon(aes(x = date, ymin = 1-Africa-Europe-Russia-NorthernAmerica-SouthernAmerica-Oceania-WesternAsia, 
                  ymax = 1-Africa-Europe-Russia-NorthernAmerica-SouthernAmerica-Oceania, fill = "7.WesternAsia"))+
  geom_ribbon(aes(x = date, ymin = 1-Africa-Europe-Russia-NorthernAmerica-SouthernAmerica-Oceania-WesternAsia-SouthernAsia, 
                  ymax = 1-Africa-Europe-Russia-NorthernAmerica-SouthernAmerica-Oceania-WesternAsia, fill = "8.SouthernAsia"))+
  geom_ribbon(aes(x = date, ymin = NorthChina+SouthChina+SoutheasternAsia, 
                  ymax = JapanKorea+NorthChina+SouthChina+SoutheasternAsia, fill = "9.Japan/Korea"))+
  geom_ribbon(aes(x = date, ymin = SouthChina+SoutheasternAsia, 
                  ymax = NorthChina+SouthChina+SoutheasternAsia, fill = "91.NorthChina"))+
  geom_ribbon(aes(x = date, ymin = SoutheasternAsia, 
                  ymax = SouthChina+SoutheasternAsia, fill = "92.SouthChina"))+
  geom_ribbon(aes(x = date, ymin = 0, 
                  ymax = SoutheasternAsia, fill = "93.SoutheasternAsia"))+
  geom_vline(xintercept = c(as.Date("2020-02-01"),as.Date("2021-08-01")), linetype = 2, linewidth = 0.2)+
  geom_line(data = trunk_even1[trunk_even1$type == "BV",],
            aes(x = date, y = Shannon_index/2.5), linetype = 4, color = "grey20")+
  scale_fill_manual("",values = value, label = label_region)+
  scale_x_date("Date", date_breaks = "6 month",date_labels = "%b %Y",
               limits = c(as.Date("2018-02-01"),as.Date("2022-08-01")), expand = c(0.01,0))+
  scale_y_continuous(expand = c(0.02,0),sec.axis = sec_axis(trans=~.* 2.5,name="Shannon index"))+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        panel.grid = element_blank())+
  # guides(fill = guide_legend(ncol = 6, nrow = 2))+
  labs(subtitle = "B/Victoria", x = "",y = "Trunk probablity")-> p3


trunk_even1_by <- dcast(trunk_even[trunk_even$type == "BY",], time ~ statistic, value.var = "mean")
trunk_even1_by$date <- decimal2Date(trunk_even1_by$time)
names(trunk_even1_by)[10] <- "JapanKorea"
ggplot(trunk_even1_by)+
  geom_ribbon(aes(x = date, ymin = 0, ymax = 1, fill = "1.Africa"))+
  geom_ribbon(aes(x = date, ymin = 1-Africa-Europe, ymax = 1-Africa, fill = "2.Europe"))+
  geom_ribbon(aes(x = date, ymin = 1-Africa-Europe-Russia, ymax = 1-Africa-Europe, fill = "3.Russia"))+
  geom_ribbon(aes(x = date, ymin = 1-Africa-Europe-Russia-NorthernAmerica, 
                  ymax = 1-Africa-Europe-Russia, fill = "4.NorthernAmerica"))+
  geom_ribbon(aes(x = date, ymin = 1-Africa-Europe-Russia-NorthernAmerica-SouthernAmerica, 
                  ymax = 1-Africa-Europe-Russia-NorthernAmerica, fill = "5.SouthernAmerica"))+
  geom_ribbon(aes(x = date, ymin = 1-Africa-Europe-Russia-NorthernAmerica-SouthernAmerica-Oceania, 
                  ymax = 1-Africa-Europe-Russia-NorthernAmerica-SouthernAmerica, fill = "6.Oceania"))+
  geom_ribbon(aes(x = date, ymin = 1-Africa-Europe-Russia-NorthernAmerica-SouthernAmerica-Oceania-WesternAsia, 
                  ymax = 1-Africa-Europe-Russia-NorthernAmerica-SouthernAmerica-Oceania, fill = "7.WesternAsia"))+
  geom_ribbon(aes(x = date, ymin = 1-Africa-Europe-Russia-NorthernAmerica-SouthernAmerica-Oceania-WesternAsia-SouthernAsia, 
                  ymax = 1-Africa-Europe-Russia-NorthernAmerica-SouthernAmerica-Oceania-WesternAsia, fill = "8.SouthernAsia"))+
  geom_ribbon(aes(x = date, ymin = NorthChina+SouthChina+SoutheasternAsia, 
                  ymax = JapanKorea+NorthChina+SouthChina+SoutheasternAsia, fill = "9.Japan/Korea"))+
  geom_ribbon(aes(x = date, ymin = SouthChina+SoutheasternAsia, 
                  ymax = NorthChina+SouthChina+SoutheasternAsia, fill = "91.NorthChina"))+
  geom_ribbon(aes(x = date, ymin = SoutheasternAsia, 
                  ymax = SouthChina+SoutheasternAsia, fill = "92.SouthChina"))+
  geom_ribbon(aes(x = date, ymin = 0, 
                  ymax = SoutheasternAsia, fill = "93.SoutheasternAsia"))+
  # geom_vline(xintercept = c(as.Date("2020-02-01"),as.Date("2021-08-01")), linetype = 2, linewidth = 0.2)+
  geom_line(data = trunk_even1[trunk_even1$type == "BY",],
            aes(x = date, y = Shannon_index/2.5), linetype = 4, color = "grey20")+
  scale_fill_manual("",values = value, label = label_region)+
  scale_x_date("Date", date_breaks = "1 year",date_labels = "%b %Y",
               limits = c(as.Date("2014-01-01"),as.Date("2019-01-01")), expand = c(0.01,0))+
  scale_y_continuous(expand = c(0.02,0),sec.axis = sec_axis(trans=~.* 2.5,name="Shannon index"))+
  theme_bw()+
  annotate("text", x = c(as.Date("2016-06-20")), y = 1.1, label = c("Before the COVID-19 pandemic"))+
  annotate("segment", x = c(as.Date("2014-01-01"),as.Date("2019-01-01")),xend = c(as.Date("2014-01-01"),as.Date("2019-01-01")),
           yend = 1.02, y = 1.17, arrow = arrow(type = "open", length = unit(0.05, "inches")))+
  coord_cartesian(ylim = c(0,1), clip = "off")+
  theme(axis.title.x = element_blank(),
        panel.border = element_rect(linewidth = 0.3),
        panel.grid = element_blank())+
  # guides(fill = guide_legend(ncol = 6, nrow = 2))+
  labs(subtitle = "B/Yamagata\n", x = "",y = "Trunk probablity")-> p4

#==output==
svg("output/Fig3.svg",width = 9, height = 8)
((p1 + p2 + p3 + p4) +plot_annotation(tag_level = "a") + plot_layout(nrow = 4,ncol = 1,guides = "collect",heights = c(1,1,1))&theme(legend.position = "right"))
dev.off()

pdf("output/Fig3.pdf",width = 9, height = 8)
((p1 + p2 + p3 + p4) +plot_annotation(tag_level = "a") + plot_layout(nrow = 4,ncol = 1,guides = "collect",heights = c(1,1,1))&theme(legend.position = "right"))
dev.off()

tiff("output/Fig3.tif",width = 9, height = 8, units = "in", res = 400, compression = "lzw")
((p1 + p2 + p3 + p4) +plot_annotation(tag_level = "a") + plot_layout(nrow = 4,ncol = 1,guides = "collect",heights = c(1,1,1))&theme(legend.position = "right"))
dev.off()

# 
# svg("output/Trunk_prop.svg",width = 9, height = 8)
# ((p1 + p2 + p3 + p4) +plot_annotation(tag_level = "a") + plot_layout(nrow = 4,ncol = 1,guides = "collect",heights = c(1,1,1))&theme(legend.position = "right"))
# dev.off()
# 
# pdf("output/Trunk_prop.pdf",width = 9, height = 8)
# ((p1 + p2 + p3 + p4) +plot_annotation(tag_level = "a") + plot_layout(nrow = 4,ncol = 1,guides = "collect",heights = c(1,1,1))&theme(legend.position = "right"))
# dev.off()
# 
# tiff("output/Trunk_prop.tif",width = 9, height = 8, units = "in", res = 400, compression = "lzw")
# ((p1 + p2 + p3 + p4) +plot_annotation(tag_level = "a") + plot_layout(nrow = 4,ncol = 1,guides = "collect",heights = c(1,1,1))&theme(legend.position = "right"))
# dev.off()
# 
# 
# svg("output/Trunk_prop_tempo.svg",width = 9, height = 8)
# ((p1 + p2 + p3 + p4) +plot_annotation(tag_level = "a") + plot_layout(nrow = 4,ncol = 1,guides = "collect",heights = c(1,1,1))&theme(legend.position = "right"))
# dev.off()
# 
# pdf("output/Trunk_prop_tempo.pdf",width = 9, height = 8)
# ((p1 + p2 + p3 + p4) +plot_annotation(tag_level = "a") + plot_layout(nrow = 4,ncol = 1,guides = "collect",heights = c(1,1,1))&theme(legend.position = "right"))
# dev.off()
# 
# tiff("output/Trunk_prop_tempo.tif",width = 9, height = 8, units = "in", res = 400, compression = "lzw")
# ((p1 + p2 + p3 + p4) +plot_annotation(tag_level = "a") + plot_layout(nrow = 4,ncol = 1,guides = "collect",heights = c(1,1,1))&theme(legend.position = "right"))
# dev.off()
# 
# svg("output/Trunk_prop_9month.svg",width = 9, height = 8)
# ((p1 + p2 + p3 + p4) +plot_annotation(tag_level = "a") + plot_layout(nrow = 4,ncol = 1,guides = "collect",heights = c(1,1,1))&theme(legend.position = "right"))
# dev.off()
# 
# pdf("output/Trunk_prop_9month.pdf",width = 9, height = 8)
# ((p1 + p2 + p3 + p4) +plot_annotation(tag_level = "a") + plot_layout(nrow = 4,ncol = 1,guides = "collect",heights = c(1,1,1))&theme(legend.position = "right"))
# dev.off()
# 
# tiff("output/Trunk_prop_9month.tif",width = 9, height = 8, units = "in", res = 400, compression = "lzw")
# ((p1 + p2 + p3 + p4) +plot_annotation(tag_level = "a") + plot_layout(nrow = 4,ncol = 1,guides = "collect",heights = c(1,1,1))&theme(legend.position = "right"))
# dev.off()
# 
# 
# svg("output/Trunk_prop_tempo_9month.svg",width = 9, height = 8)
# ((p1 + p2 + p3 + p4) +plot_annotation(tag_level = "a") + plot_layout(nrow = 4,ncol = 1,guides = "collect",heights = c(1,1,1))&theme(legend.position = "right"))
# dev.off()
# 
# pdf("output/Trunk_prop_tempo_9month.pdf",width = 9, height = 8)
# ((p1 + p2 + p3 + p4) +plot_annotation(tag_level = "a") + plot_layout(nrow = 4,ncol = 1,guides = "collect",heights = c(1,1,1))&theme(legend.position = "right"))
# dev.off()
# 
# tiff("output/Trunk_prop_tempo_9month.tif",width = 9, height = 8, units = "in", res = 400, compression = "lzw")
# ((p1 + p2 + p3 + p4) +plot_annotation(tag_level = "a") + plot_layout(nrow = 4,ncol = 1,guides = "collect",heights = c(1,1,1))&theme(legend.position = "right"))
# dev.off()