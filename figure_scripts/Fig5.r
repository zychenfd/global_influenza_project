#==load package==
library(ggridges)
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
library(lubridate)
library(circlize)
library(coda)
library(cowplot)
library(scales)
library(zoo)
library(grid)
library(patchwork)

#==define color==
colors <- c(pal_npg("nrc", alpha =1)(10)[c(1:7,9:10)],"darkred","#FADDA9","grey80")
colors1 <- c(pal_aaas("default", alpha =0.7)(10))
show_col(colors)
show_col(colors1)
value = c("Japan/Korea" = colors[1],"Western Asia" = colors[3],"WesternAsia" = colors[3],
          "Northern America" = colors[6],"North Am" = colors[6],"NorthernAmerica" = colors[6],
          "South-eastern Asia"= colors[4], "Southern Asia"= colors[5],"SoutheasternAsia"= colors[4], "SouthernAsia"= colors[5],
          "Europe"= colors[2], "Oceania"= colors[7],"NorthChina" = colors[8], "SouthChina" = colors[11],
          "North China" = colors[8], "South China" = colors[11], "China (N)" = colors[8], "China (S)" = colors[11],
          "Russia"= colors[10],  "Southern America"= colors[12],  "South Am"= colors[12],  "SouthernAmerica"= colors[12],
          "Africa"= colors[9], "Americas" = colors1[1],"Asia" = "#D8BFD8", "China" = colors1[4])

#==tip associated persistence==
per_h1n1_11r_e1 <- read.delim("../genomic_part/post-analyses/persistence/tip_persistence/h1n1_persistence_11region/epoch1/out.stats") %>%
  mutate(type = "H1N1", epoch = "epoch1", geo_num = "11_region")
per_h1n1_11r_e2 <- read.delim("../genomic_part/post-analyses/persistence/tip_persistence/h1n1_persistence_11region/epoch2/out.stats") %>%
  mutate(type = "H1N1", epoch = "epoch2", geo_num = "11_region")
per_h1n1_11r_e3 <- read.delim("../genomic_part/post-analyses/persistence/tip_persistence/h1n1_persistence_11region/epoch3/out.stats") %>%
  mutate(type = "H1N1", epoch = "epoch3", geo_num = "11_region")

per_h3n2_11r_e1 <- read.delim("../genomic_part/post-analyses/persistence/tip_persistence/h3n2_persistence_11region/epoch1/out.stats") %>%
  mutate(type = "H3N2", epoch = "epoch1", geo_num = "11_region")
per_h3n2_11r_e2 <- read.delim("../genomic_part/post-analyses/persistence/tip_persistence/h3n2_persistence_11region/epoch2/out.stats") %>%
  mutate(type = "H3N2", epoch = "epoch2", geo_num = "11_region")
per_h3n2_11r_e3 <- read.delim("../genomic_part/post-analyses/persistence/tip_persistence/h3n2_persistence_11region/epoch3/out.stats") %>%
  mutate(type = "H3N2", epoch = "epoch3", geo_num = "11_region")

per_bv_11r_e1 <- read.delim("../genomic_part/post-analyses/persistence/tip_persistence/bv_persistence_11region/epoch1/out.stats") %>%
  mutate(type = "BV", epoch = "epoch1", geo_num = "11_region")
per_bv_11r_e2 <- read.delim("../genomic_part/post-analyses/persistence/tip_persistence/bv_persistence_11region/epoch2/out.stats") %>%
  mutate(type = "BV", epoch = "epoch2", geo_num = "11_region")
per_bv_11r_e3 <- read.delim("../genomic_part/post-analyses/persistence/tip_persistence/bv_persistence_11region/epoch3/out.stats") %>%
  mutate(type = "BV", epoch = "epoch3", geo_num = "11_region")

per_by_11r_e1 <- read.delim("../genomic_part/post-analyses/persistence/tip_persistence/by_persistence_11region/epoch1/out.stats") %>%
  mutate(type = "BY", epoch = "epoch1", geo_num = "11_region")

per_11r <- rbind(per_h1n1_11r_e1,per_h1n1_11r_e2,per_h1n1_11r_e3,
                 per_h3n2_11r_e1,per_h3n2_11r_e2,per_h3n2_11r_e3,
                 per_bv_11r_e1,per_bv_11r_e2,per_bv_11r_e3,per_by_11r_e1)
remove(per_h1n1_11r_e1,per_h1n1_11r_e2,per_h1n1_11r_e3,
       per_h3n2_11r_e1,per_h3n2_11r_e2,per_h3n2_11r_e3,
       per_bv_11r_e1,per_bv_11r_e2,per_bv_11r_e3,per_by_11r_e1)

per_11r$statistic <- str_remove_all(per_11r$statistic, "persistence_")
factor(per_11r$statistic, levels = unique(per_11r$statistic)[rev(c(1,6,10,4,8,2,3,5,9,11,12,7))]) -> per_11r$statistic
ggplot(data = per_11r[per_11r$type == "H1N1",])+
  geom_hline(yintercept = 0.5, color = "red", linetype = 2, size = 0.3)+
  geom_errorbar(aes(x = statistic, ymin = lower, ymax = upper, color = epoch),position = position_dodge(width = 0.5),width = 0.25)+
  geom_point(aes(x = statistic, y = mean, fill = epoch), shape =21, position = position_dodge(width = 0.5))+
  scale_y_continuous(limits = c(0,2.5))+
  theme_bw()+
  coord_flip()+
  theme(legend.position = c(0.7,0.8),
        panel.grid.minor = element_blank())+
  scale_fill_manual("Epochs",values = colors1, labels = c("Epoch 1","Epoch 2","Epoch 3"))+
  scale_color_manual("Epochs",values = colors1, labels = c("Epoch 1","Epoch 2","Epoch 3"))+
  scale_x_discrete(labels = rev(c("Global", "North America", "South America","Europe" , "Russia",
                              "Africa", "China","Japan/Korea","Southeast Asia",  "South Asia",
                              "West Asia", "Oceania" )))+
  theme(axis.title.y = element_blank())+
  labs(subtitle = "H1N1pdm09", y = "Tips-related persistence (years)") -> p1

ggplot(data = per_11r[per_11r$type == "H3N2",])+
  geom_hline(yintercept = 0.5, color = "red", linetype = 2, size = 0.3)+
  geom_errorbar(aes(x = statistic, ymin = lower, ymax = upper, color = epoch),position = position_dodge(width = 0.5),width = 0.25)+
  geom_point(aes(x = statistic, y = mean, fill = epoch), shape =21, position = position_dodge(width = 0.5))+
  scale_y_continuous(limits = c(0,2.5))+
  theme_bw()+
  coord_flip()+
  guides(fill = F, color = F)+
  scale_fill_manual(values = colors1)+
  scale_color_manual(values = colors1)+
  theme(axis.title.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank())+
  labs(subtitle = "H3N2", y = "Tips-related persistence (years)") -> p2

ggplot(data = per_11r[per_11r$type == "BV",])+
  geom_hline(yintercept = 0.5, color = "red", linetype = 2, size = 0.3)+
  geom_errorbar(aes(x = statistic, ymin = lower, ymax = upper, color = epoch),position = position_dodge(width = 0.5),width = 0.25)+
  geom_point(aes(x = statistic, y = mean, fill = epoch), shape =21, position = position_dodge(width = 0.5))+
  scale_y_continuous(limits = c(0,2.5))+
  theme_bw()+
  coord_flip()+
  guides(fill = F, color = F)+
  scale_fill_manual(values = colors1)+
  scale_color_manual(values = colors1)+
  theme(axis.title.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank())+
  labs(subtitle = "B/Victoria", y = "Tips-related persistence (years)")-> p3

ggplot(data = per_11r[per_11r$type == "BY",])+
  geom_hline(yintercept = 0.5, color = "red", linetype = 2, size = 0.3)+
  geom_errorbar(aes(x = statistic, ymin = lower, ymax = upper),width = 0.125, color = colors1[1])+
  geom_point(aes(x = statistic, y = mean), shape =21, fill = colors1[1])+
  scale_y_continuous(limits = c(0,2.5))+
  coord_flip()+
  theme_bw()+
  theme(axis.title.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank())+
  labs(subtitle = "B/Yamagata", y = "Tips-related persistence (years)") -> p4

#===lineage associated persistence (H3N2)==
h3n2_p <- read.csv("../genomic_part/post-analyses/persistence/lineage_persistence/h3n2_persistence.csv")
h3n2_p$evaluationDate <- as.Date(date_decimal(decimal_date(as.Date("2023-07-27")) - h3n2_p$evaluationTime))
h3n2_p2 <- h3n2_p[order(h3n2_p$stateAtEvaluationTime, h3n2_p$evaluationDate),]
h3n2_p3 <- h3n2_p2 %>% 
  group_by(stateAtEvaluationTime) %>%
  mutate(Persistent_mean_roll = rollmean(Persistent_mean, k=3, fill=NA, na.pad=T, align='center'))

ggplot(data = h3n2_p3[h3n2_p3$stateAtEvaluationTime %in% c("SoutheasternAsia","Africa","SouthernAsia") &
                        h3n2_p3$evaluationDate >= as.Date("2017-01-01") & h3n2_p3$evaluationDate <= as.Date("2023-03-31"),]) + 
  annotate("rect", xmin = as.Date("2020-02-01"),xmax = as.Date("2021-07-31"),
           ymin = 0,ymax = 2.4,alpha = 0.2,fill = colors[5])+
  geom_line(aes(x = evaluationDate, y = Persistent_mean_roll, color = stateAtEvaluationTime))+
  geom_point(aes(x = evaluationDate, y = Persistent_mean_roll, fill = stateAtEvaluationTime, size = No_lineage), shape = 21, stroke = 0.1)+
  scale_x_date("Evaluation date",date_labels = "%b %Y", date_breaks = "6 months")+
  scale_y_continuous("Lineage-related persistence (years)", limits = c(-0.01,2.4),expand = c(0,0))+
  scale_color_manual("Geographic regions", 
                     values = value,
                     labels = c("Africa","Southeast Asia","South Asia"))+
  scale_fill_manual("Geographic regions",
                    values = value,
                    labels = c("Africa", "Southeast Asia","South Asia"))+
  theme_bw()+
  scale_size_continuous("Number of lineages at each timepoint")+
  theme(panel.grid = element_blank(),
        legend.position = c(0.2,0.65),
        panel.background = element_rect(fill = "transparent", color = "transparent"),
        legend.background = element_blank(),
        plot.background = element_rect(fill = "transparent", color = "transparent"))+
  guides(fill = guide_colorbar(frame.colour = "black",
                               title.position="top",
                               ticks.colour = "black"),
         color = guide_legend(title.position="top", nrow = 1),
         fill = guide_legend(title.position="top", nrow = 1),
         size = guide_legend(title.position="top", nrow = 1))+
  labs(subtitle = "H3N2", tag = "e") -> p5

#==plot for result of regression model== 
ag <- read.csv("../model_part/antigenic_effect_sumary.csv") %>% filter(year != 2017)
npi <- read.csv("../model_part/air-traffic_effect_sumary.csv") %>% filter(year != 2017)

factor(ag$year, levels = unique(ag$year)) -> ag$year
factor(npi$year, levels = unique(npi$year)) -> npi$year
ggplot(data = ag[ag$region == "Africa",])+
  geom_hline(yintercept = 0, color = "red", linetype = 2, size = 0.3)+
  geom_errorbar(aes(x = year, ymin = hdi_5., ymax = hdi_95.), width = 0.1, color = colors[9])+
  geom_point(aes(x = year, y = mean), shape = 21, fill = colors[9], size = 2)+
  coord_flip(ylim = c(-0.8,1.5))+
  theme_bw()+
  theme(plot.subtitle = element_text(hjust = 0.5),
        panel.grid.minor = element_blank(),
        plot.margin =  margin(0, 0, 0, 0, "cm"))+
  labs(x = "Years", y = "Effect of antigenic drift on persistence",
       subtitle = "Africa", tag = "f") -> p6

ggplot(data = npi[npi$region == "Africa",])+
  geom_hline(yintercept = 0, color = "red", linetype = 2, size = 0.3)+
  geom_errorbar(aes(x = year, ymin = hdi_5., ymax = hdi_95.), width = 0.1, color = colors[9])+
  geom_point(aes(x = year, y = mean), shape = 21, fill = colors[9], size = 2)+
  coord_flip(ylim = c(-0.8,1.5))+
  theme_bw()+
  theme(plot.margin =  margin(0.2, 0, 0, 0, "cm"),
        panel.grid.minor = element_blank(),)+
  labs(x = "Years", y = "Effect of air traffic on persistence") -> p7

ggplot(data = ag[ag$region == "Southern Asia",])+
  geom_hline(yintercept = 0, color = "red", linetype = 2, size = 0.3)+
  geom_errorbar(aes(x = year, ymin = hdi_5., ymax = hdi_95.), width = 0.1, color = colors[5])+
  geom_point(aes(x = year, y = mean), shape = 21, fill = colors[5], size = 2)+
  coord_flip(ylim = c(-1.2,2))+
  theme_bw()+
  theme(plot.subtitle = element_text(hjust = 0.5),
        panel.grid.minor = element_blank(),
        plot.margin =  margin(0, 0, 0, 0, "cm"),
        axis.text.y = element_blank(),
        # axis.title.y = element_blank()
        )+
  labs(x = "", y = "Effect of antigenic drift on persistence",
       subtitle = "South Asia") -> p8

ggplot(data = npi[npi$region == "Southern Asia",])+
  geom_hline(yintercept = 0, color = "red", linetype = 2, size = 0.3)+
  geom_errorbar(aes(x = year, ymin = hdi_5., ymax = hdi_95.), width = 0.1, color = colors[5])+
  geom_point(aes(x = year, y = mean), shape = 21, fill = colors[5], size = 2)+
  coord_flip(ylim = c(-1.2,2))+
  theme_bw()+
  theme( 
    # axis.title.y = element_blank(),
         panel.grid.minor = element_blank(),
         axis.text.y = element_blank(),
         plot.margin =  margin(0.2, 0, 0, 0, "cm"))+
  labs(x = "", y = "Effect of air traffic on persistence") -> p9

ggplot(data = ag[ag$region == "South-eastern Asia",])+
  geom_hline(yintercept = 0, color = "red", linetype = 2, size = 0.3)+
  geom_errorbar(aes(x = year, ymin = hdi_5., ymax = hdi_95.), width = 0.1, color = colors[4])+
  geom_point(aes(x = year, y = mean), shape = 21, fill = colors[4], size = 2)+
  coord_flip(ylim = c(-0.8,1.5))+
  theme_bw()+
  theme(plot.subtitle = element_text(hjust = 0.5),
        plot.margin =  margin(0, 0, 0, 0, "cm"),
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        # axis.title.y = element_blank()
        )+
  labs(x = "", y = "Effect of antigenic drift on persistence",
       subtitle = "Southeast Asia") -> p10

ggplot(data = npi[npi$region == "South-eastern Asia",])+
  geom_hline(yintercept = 0, color = "red", linetype = 2, size = 0.3)+
  geom_errorbar(aes(x = year, ymin = hdi_5., ymax = hdi_95.), width = 0.1, color = colors[4])+
  geom_point(aes(x = year, y = mean), shape = 21, fill = colors[4], size = 2)+
  coord_flip(ylim = c(-0.8,1.5))+
  theme_bw()+
  theme( 
    # axis.title.y = element_blank(),
         axis.text.y = element_blank(),
         panel.grid.minor = element_blank(),
         plot.margin =  margin(0.2, 0, 0, 0, "cm") )+
  labs(x = "", y = "Effect of air traffic on persistence")-> p11

#==output==
pdf("output/Fig5.pdf",width = 10.5, height = 10)
(p1|p2|p3|p4)+plot_annotation(tag_levels = "a") -> f1
p5 -> f2
((p6/p7)|(p8/p9)|(p10/p11)) -> f3
viewport(x = 0, y = 0.56, width = 1, height = 0.45, just = c("left", "bottom")) -> vp1
viewport(x = 0, y = 0.29, width = 1, height = 0.3, just = c("left", "bottom")) -> vp2
viewport(x = 0, y = -0.005, width = 1, height = 0.305, just = c("left", "bottom")) -> vp3
print(f1,vp = vp1)
print(f2, vp = vp2)
print(f3, vp = vp3)
dev.off()

svg("output/Fig5.svg",width = 10.5, height = 10)
(p1|p2|p3|p4)+plot_annotation(tag_levels = "a") -> f1
p5 -> f2
((p6/p7)|(p8/p9)|(p10/p11)) -> f3
viewport(x = 0, y = 0.56, width = 1, height = 0.45, just = c("left", "bottom")) -> vp1
viewport(x = 0, y = 0.29, width = 1, height = 0.3, just = c("left", "bottom")) -> vp2
viewport(x = 0, y = -0.005, width = 1, height = 0.305, just = c("left", "bottom")) -> vp3
print(f1,vp = vp1)
print(f2, vp = vp2)
print(f3, vp = vp3)
dev.off()