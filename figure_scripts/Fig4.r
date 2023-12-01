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
library(lubridate)
library(circlize)
library(coda)
library(cowplot)
library(scales)
library(networkD3)
library(grid)

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
          "Africa"= colors[9], "Americas" = colors1[1], "Asia" = "#D8BFD8", "China" = colors1[4])

#=============
#==h1n1pdm09==
#=============
h1n1_jump <- read.delim("../genomic_part/post-analyses/jump_history/h1n1_even_Jumps.txt", sep = ",")
h1n1_jump$date <- decimal2Date(decimal_date(as.Date("2023-07-31")) - h1n1_jump$time)
h1n1_jump$period[h1n1_jump$date < as.Date("2020-02-01") & h1n1_jump$date >= as.Date("2018-02-01")] <- "Period 1"
h1n1_jump$period[h1n1_jump$date < as.Date("2021-08-01") & h1n1_jump$date >= as.Date("2020-02-01")] <- "Period 2"
h1n1_jump$period[h1n1_jump$date >= as.Date("2021-08-01")] <- "Period 3"
h1n1_jump1 <- h1n1_jump %>%
  filter(!is.na(period)) %>%
  group_by(treeId, period, startLocation, endLocation) %>%
  summarise(n = n())
h1n1_jump2 <- h1n1_jump1 %>%
  group_by(period, startLocation, endLocation) %>%
  summarise(mean = mean(n),
            sd = sd(n))
h1n1_jump2$startLocation[h1n1_jump2$startLocation == "SouthChina"] <- "China (S)"
h1n1_jump2$startLocation[h1n1_jump2$startLocation == "SoutheasternAsia"] <- "South-eastern Asia"
h1n1_jump2$startLocation[h1n1_jump2$startLocation == "SouthernAmerica"] <- "South Am"
h1n1_jump2$startLocation[h1n1_jump2$startLocation == "SouthernAsia"] <- "Southern Asia"
h1n1_jump2$startLocation[h1n1_jump2$startLocation == "WesternAsia"] <- "Western Asia"
h1n1_jump2$startLocation[h1n1_jump2$startLocation == "NorthChina"] <- "China (N)"
h1n1_jump2$startLocation[h1n1_jump2$startLocation == "NorthernAmerica"] <- "North Am"
h1n1_jump2$endLocation[h1n1_jump2$endLocation == "SouthChina"] <- "China (S)"
h1n1_jump2$endLocation[h1n1_jump2$endLocation == "SoutheasternAsia"] <- "South-eastern Asia"
h1n1_jump2$endLocation[h1n1_jump2$endLocation == "SouthernAmerica"] <- "South Am"
h1n1_jump2$endLocation[h1n1_jump2$endLocation == "SouthernAsia"] <- "Southern Asia"
h1n1_jump2$endLocation[h1n1_jump2$endLocation == "WesternAsia"] <- "Western Asia"
h1n1_jump2$endLocation[h1n1_jump2$endLocation == "NorthChina"] <- "China (N)"
h1n1_jump2$endLocation[h1n1_jump2$endLocation == "NorthernAmerica"] <- "North Am"
levels <- c("South-eastern Asia","China (N)","China (S)","Japan/Korea", "Southern Asia", "Western Asia",
            "North Am", "South Am",  "Europe",  "Russia", "Africa",  "Oceania")

sum(h1n1_jump2$mean[h1n1_jump2$period == "Period 1"])
sum(h1n1_jump2$mean[h1n1_jump2$period == "Period 2"])
sum(h1n1_jump2$mean[h1n1_jump2$period == "Period 3"])

#sankey plot (epoch 1)
data_tmp <- h1n1_jump2[h1n1_jump2$period == "Period 1",] %>% group_by(startLocation) %>% summarise(n = sum(mean))
data_tmp1 <- h1n1_jump2[h1n1_jump2$period == "Period 1",] %>% group_by(endLocation) %>% summarise(n = sum(mean))
name <- data.frame(name = c(data_tmp$startLocation[order(data_tmp$n, decreasing = T)],
                            data_tmp1$endLocation[order(data_tmp1$n, decreasing = T)]))
names1 <- data.frame(location = name$name[1:12], source = 0:11)
names2 <- data.frame(location = name$name[13:24], target = 12:23)
h1n1_jump3 <- left_join(h1n1_jump2[h1n1_jump2$period == "Period 1",], names1, by = c("startLocation" = "location"))
h1n1_jump3 <- left_join(h1n1_jump3, names2, by = c("endLocation" = "location"))
name$name <- str_remove_all(name$name," |/|-")
my_color <- 'd3.scaleOrdinal() .domain(["China(N)","SoutheasternAsia","China(S)","JapanKorea", "SouthernAsia", "WesternAsia","NorthAm", "SouthAm",  "Europe",  "Russia", "Africa",  "Oceania"]) .range(["#7E6148FF","#3C5488FF","#FADDA9","#E64B35FF", "#F39B7FFF" , "#00A087FF" ,"#8491B4FF","#CCCCCC","#4DBBD5FF","darkred","#B09C85FF" ,"#91D1C2FF"])'

sankeyNetwork(Links = data.frame(h1n1_jump3), Nodes = name,iterations = 0,
              colourScale = my_color,
              Source = "source", Target = "target", Value = "mean", NodeID = "name",fontFamily = NULL,
              sinksRight=FALSE, nodeWidth=10, nodePadding=5,height = 400, fontSize = 0) 

#sankey plot (epoch 2)
data_tmp <- h1n1_jump2[h1n1_jump2$period == "Period 2",] %>% group_by(startLocation) %>% summarise(n = sum(mean))
data_tmp1 <- h1n1_jump2[h1n1_jump2$period == "Period 2",] %>% group_by(endLocation) %>% summarise(n = sum(mean))
name_epoch2 <- data.frame(name = c(data_tmp$startLocation[order(data_tmp$n, decreasing = T)],
                                   data_tmp1$endLocation[order(data_tmp1$n, decreasing = T)]))
names1 <- data.frame(location = name_epoch2$name[1:12], source = 0:11)
names2 <- data.frame(location = name_epoch2$name[13:24], target = 12:23)
h1n1_jump4 <- left_join(h1n1_jump2[h1n1_jump2$period == "Period 2",], names1, by = c("startLocation" = "location"))
h1n1_jump4 <- left_join(h1n1_jump4, names2, by = c("endLocation" = "location"))
name_epoch2$name <- str_remove_all(name_epoch2$name," |/|-")

sankeyNetwork(Links = data.frame(h1n1_jump4), Nodes = name_epoch2,iterations = 0,
              colourScale = my_color,
              Source = "source", Target = "target", Value = "mean", NodeID = "name",fontFamily = NULL,
              sinksRight=FALSE, nodeWidth=10, nodePadding=5,height = 400, fontSize = 0)

#sankey plot (epoch 3)
data_tmp <- h1n1_jump2[h1n1_jump2$period == "Period 3",] %>% group_by(startLocation) %>% summarise(n = sum(mean))
data_tmp1 <- h1n1_jump2[h1n1_jump2$period == "Period 3",] %>% group_by(endLocation) %>% summarise(n = sum(mean))
name_epoch3 <- data.frame(name = c(data_tmp$startLocation[order(data_tmp$n, decreasing = T)],
                            data_tmp1$endLocation[order(data_tmp1$n, decreasing = T)]))
names1 <- data.frame(location = name_epoch3$name[1:12], source = 0:11)
names2 <- data.frame(location = name_epoch3$name[13:24], target = 12:23)
h1n1_jump5 <- left_join(h1n1_jump2[h1n1_jump2$period == "Period 3",], names1, by = c("startLocation" = "location"))
h1n1_jump5 <- left_join(h1n1_jump5, names2, by = c("endLocation" = "location"))
name_epoch3$name <- str_remove_all(name_epoch3$name," |/|-")

sankeyNetwork(Links = data.frame(h1n1_jump5), Nodes = name_epoch3,iterations = 0,
              colourScale = my_color,
              Source = "source", Target = "target", Value = "mean", NodeID = "name",fontFamily = NULL,
              sinksRight=FALSE, nodeWidth=10, nodePadding=5,height = 400, fontSize = 0) 

#==Net export plot==
h1n1_jump3 <- h1n1_jump1 %>%
  group_by(treeId, period, startLocation) %>%
  summarise(export_n = sum(n))
h1n1_jump3_1 <- h1n1_jump1 %>%
  group_by(treeId, period, endLocation) %>%
  summarise(import_n = sum(n))
h1n1_jump4 <- full_join(h1n1_jump3, h1n1_jump3_1, by = c("treeId" = "treeId", "period" = "period", "startLocation" = "endLocation"))
h1n1_jump4[is.na(h1n1_jump4)] <- 0
h1n1_jump5 <- h1n1_jump4 %>%
  mutate(net_export = export_n - import_n) %>%
  group_by(period, startLocation) %>%
  summarise(net_export_mean = mean(net_export),
            net_export_sd = sd(net_export),
            net_export_1q = quantile(net_export,0.25),
            net_export_3q = quantile(net_export,0.75))

h1n1_jump5 <- h1n1_jump5[order(h1n1_jump5$net_export_mean, decreasing = F),]
h1n1_jump5_1 <- h1n1_jump5[h1n1_jump5$period == "Period 1",]
factor(h1n1_jump5_1$startLocation, levels = h1n1_jump5_1$startLocation) -> h1n1_jump5_1$startLocation
ggplot(data = h1n1_jump5_1) +
  geom_hline(yintercept = 0, linetype = 2)+
  geom_bar(aes(x = startLocation, y = net_export_mean, fill = startLocation), stat = "identity", width = 0.5)+
  geom_errorbar(aes(x = startLocation, ymin = net_export_1q, ymax = net_export_3q), width = 0.2, size = 0.2)+
  theme_bw()+
  theme(axis.text.y = element_blank(),
        plot.margin =  margin(0.1, 0.1, 0.1, 0.1, "cm"))+
  coord_flip()+
  scale_y_continuous(limits = c(min(h1n1_jump5$net_export_1q), max(h1n1_jump5$net_export_3q)), breaks = seq(-30,30,30))+
  scale_fill_manual(values = value)+
  labs(x = "", y = "Net exports",subtitle = "") +
  guides(fill = F)-> fig1_1

h1n1_jump5_2 <- h1n1_jump5[h1n1_jump5$period == "Period 2",]
factor(h1n1_jump5_2$startLocation, levels = h1n1_jump5_2$startLocation) -> h1n1_jump5_2$startLocation
ggplot(data = h1n1_jump5_2) +
  geom_hline(yintercept = 0, linetype = 2)+
  geom_bar(aes(x = startLocation, y = net_export_mean, fill = startLocation), stat = "identity", width = 0.5)+
  geom_errorbar(aes(x = startLocation, ymin = net_export_1q, ymax = net_export_3q), width = 0.2, size = 0.2)+
  theme_bw()+
  theme(axis.text.y = element_blank(),
        plot.margin =  margin(0.1, 0.1, 0.1, 0.1, "cm"))+
  coord_flip()+
  scale_y_continuous(limits = c(min(h1n1_jump5$net_export_1q), max(h1n1_jump5$net_export_3q)), breaks = seq(-30,30,30))+
  scale_fill_manual(values = c( "Japan/Korea" = colors[4],"WesternAsia" = colors[6],"NorthernAmerica" = colors[7],
                                "SoutheasternAsia"= colors[1], "SouthernAsia"= colors[5], "Europe"= colors[9], "Oceania"= colors[12],
                                "NorthChina" = colors[2], "SouthChina" = colors[3], "Russia"= colors[10],  "SouthernAmerica"= colors[8], "Africa"= colors[11]))+
  labs(x = "", y = "Net exports",subtitle = "") +
  guides(fill = F) -> fig2_1

h1n1_jump5_3 <- h1n1_jump5[h1n1_jump5$period == "Period 3",]
factor(h1n1_jump5_3$startLocation, levels = h1n1_jump5_3$startLocation) -> h1n1_jump5_3$startLocation
ggplot(data = h1n1_jump5_3) +
  geom_hline(yintercept = 0, linetype = 2)+
  geom_bar(aes(x = startLocation, y = net_export_mean, fill = startLocation), stat = "identity", width = 0.5)+
  geom_errorbar(aes(x = startLocation, ymin = net_export_1q, ymax = net_export_3q), width = 0.2, size = 0.2)+
  theme_bw()+
  theme(axis.text.y = element_blank(),
        plot.margin =  margin(0.1, 0.1, 0.1, 0.1, "cm"))+
  coord_flip()+
  scale_y_continuous(limits = c(min(h1n1_jump5$net_export_1q), max(h1n1_jump5$net_export_3q)), breaks = seq(-30,30,30))+
  scale_fill_manual(values = c( "Japan/Korea" = colors[4],"WesternAsia" = colors[6],"NorthernAmerica" = colors[7],
                                "SoutheasternAsia"= colors[1], "SouthernAsia"= colors[5], "Europe"= colors[9], "Oceania"= colors[12],
                                "NorthChina" = colors[2], "SouthChina" = colors[3],  "Russia"= colors[10],  "SouthernAmerica"= colors[8], "Africa"= colors[11]))+
  labs(x = "", y = "Net exports",subtitle = "") +
  guides(fill = F) -> fig3_1

#legend panel
ggplot(data = h1n1_jump5_1) +
  geom_bar(aes(x = startLocation, y = 0, fill = startLocation), stat = "identity", width = 0.5)+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.background = element_rect(fill = "transparent",color = "transparent"),
        panel.background = element_rect(fill = "transparent",color = "transparent"),
        panel.border = element_rect(fill = "transparent",color = "transparent"),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12),
        plot.margin =  margin(0.1, 0.1, 0.1, 0.1, "cm"),
        panel.grid = element_blank(),
        legend.position = c(0.5,0.5))+
  scale_fill_manual("Geographic regions",
                    values = value,
                    labels = c("Africa","South America","Russia","North China","South China",
                               "South Asia","Oceania", "Europe", "North America", "West Asia",
                               "Southeast Asia","Japan/Korea"))+
  labs(x = "", y = "",title = "") +
  guides(fill = guide_legend(nrow = 3))-> fig_legend

#blank panel
ggplot(data = h1n1_jump5_1) +
  geom_bar(aes(x = startLocation, y = 0, fill = startLocation), stat = "identity", width = 0.5)+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent",color = "transparent"),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12),
        plot.margin =  margin(0.1, 0.1, 0.1, 0.1, "cm"),
        panel.grid = element_blank(),
        legend.position = "none")+
  scale_fill_manual("Geographic regions",
                    values = value,
                    labels = c("Africa","South America","Russia","North China","South China",
                               "South Asia","Oceania", "Europe", "North America", "West Asia",
                               "Southeast Asia","Japan/Korea"))+
  labs(x = "", y = "",title = "") +
  guides(fill = guide_legend(nrow = 3))-> fig_white

#panel a (overall markov jumps)
date_match <- data.frame(date_week = seq(as.Date("1995-01-02"),as.Date("2023-12-25"),7))
date_match$ISO_YEAR <- isoyear(date_match$date)
date_match$ISO_WEEK <- isoweek(date_match$date)
h1n1_jump6 <- as.data.frame(h1n1_jump) %>%
  mutate(ISO_YEAR = isoyear(date),
         ISO_WEEK = isoweek(date),
         treeId = as.character(treeId)) %>%
  filter(!is.na(period)) %>%
  group_by(treeId, ISO_YEAR, ISO_WEEK) %>%
  summarise(n = n())
h1n1_jump6 <- h1n1_jump6 %>%
  left_join(date_match) %>%
  group_by(date_week) %>%
  summarise(mean = mean(n),
            sd = sd(n))

ggplot(data = h1n1_jump6) +
  # geom_vline(xintercept = c(as.Date("2020-02-01"),as.Date("2021-08-01")), linetype = 2, size = 0.2)+
  annotate("rect", xmin = as.Date("2020-02-01"),xmax = as.Date("2021-07-31"),
           ymin = 0,ymax = 25,alpha = 0.2,fill = colors[5])+
  geom_ribbon(aes(x = date_week, ymin = 0, ymax = mean), fill = colors1[1], alpha = 0.7)+
  geom_line(aes(x = date_week, y = mean, group = 1))+
  scale_y_continuous("Markov jumps", limits = c(-0.5,25), expand = c(0,0), breaks = seq(0,25,10))+
  scale_x_date(limits = c(as.Date("2018-02-01"), as.Date("2023-07-31")), expand = c(0.01,0))+
  theme_bw()+
  annotate("text", x = as.Date("2018-03-01"), y = 20, label = "H1N1pdm09", hjust = 0, size = 4)+
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        plot.margin =  margin(0.1, 0.1, 0.1, 0.1, "cm"),
        axis.text = element_text(size = 12),
        panel.grid.minor = element_blank()) -> m1

#========
#==h3n2==
#========
h3n2_jump <- read.delim("../genomic_part/post-analyses/jump_history/h3n2_even_Jumps.txt", sep = ",")
h3n2_jump$date <- decimal2Date(decimal_date(as.Date("2023-07-27")) - h3n2_jump$time)
h3n2_jump$period[h3n2_jump$date < as.Date("2020-02-01") & h3n2_jump$date >= as.Date("2018-02-01")] <- "Period 1"
h3n2_jump$period[h3n2_jump$date < as.Date("2021-08-01") & h3n2_jump$date >= as.Date("2020-02-01")] <- "Period 2"
h3n2_jump$period[h3n2_jump$date >= as.Date("2021-08-01")] <- "Period 3"
h3n2_jump1 <- h3n2_jump %>%
  filter(!is.na(period)) %>%
  group_by(treeId, period, startLocation, endLocation) %>%
  summarise(n = n())
h3n2_jump2 <- h3n2_jump1 %>%
  group_by(period, startLocation, endLocation) %>%
  summarise(mean = mean(n),
            sd = mean(n))
h3n2_jump2$startLocation[h3n2_jump2$startLocation == "SouthChina"] <- "China (S)"
h3n2_jump2$startLocation[h3n2_jump2$startLocation == "SoutheasternAsia"] <- "South-eastern Asia"
h3n2_jump2$startLocation[h3n2_jump2$startLocation == "SouthernAmerica"] <- "South Am"
h3n2_jump2$startLocation[h3n2_jump2$startLocation == "SouthernAsia"] <- "Southern Asia"
h3n2_jump2$startLocation[h3n2_jump2$startLocation == "WesternAsia"] <- "Western Asia"
h3n2_jump2$startLocation[h3n2_jump2$startLocation == "NorthChina"] <- "China (N)"
h3n2_jump2$startLocation[h3n2_jump2$startLocation == "NorthernAmerica"] <- "North Am"
h3n2_jump2$endLocation[h3n2_jump2$endLocation == "SouthChina"] <- "China (S)"
h3n2_jump2$endLocation[h3n2_jump2$endLocation == "SoutheasternAsia"] <- "South-eastern Asia"
h3n2_jump2$endLocation[h3n2_jump2$endLocation == "SouthernAmerica"] <- "South Am"
h3n2_jump2$endLocation[h3n2_jump2$endLocation == "SouthernAsia"] <- "Southern Asia"
h3n2_jump2$endLocation[h3n2_jump2$endLocation == "WesternAsia"] <- "Western Asia"
h3n2_jump2$endLocation[h3n2_jump2$endLocation == "NorthChina"] <- "China (N)"
h3n2_jump2$endLocation[h3n2_jump2$endLocation == "NorthernAmerica"] <- "North Am"

#sankey plot (epoch 1)
data_tmp <- h3n2_jump2[h3n2_jump2$period == "Period 1",] %>% group_by(startLocation) %>% summarise(n = sum(mean))
data_tmp1 <- h3n2_jump2[h3n2_jump2$period == "Period 1",] %>% group_by(endLocation) %>% summarise(n = sum(mean))
name <- data.frame(name = c(data_tmp$startLocation[order(data_tmp$n, decreasing = T)],
                            data_tmp1$endLocation[order(data_tmp1$n, decreasing = T)]))
names1 <- data.frame(location = name$name[1:12], source = 0:11)
names2 <- data.frame(location = name$name[13:24], target = 12:23)
h3n2_jump3 <- left_join(h3n2_jump2[h3n2_jump2$period == "Period 1",], names1, by = c("startLocation" = "location"))
h3n2_jump3 <- left_join(h3n2_jump3, names2, by = c("endLocation" = "location"))
name$name <- str_remove_all(name$name," |/|-")
my_color <- 'd3.scaleOrdinal() .domain(["China(N)","SoutheasternAsia","China(S)","JapanKorea", "SouthernAsia", "WesternAsia","NorthAm", "SouthAm",  "Europe",  "Russia", "Africa",  "Oceania"]) .range(["#7E6148FF","#3C5488FF","#FADDA9","#E64B35FF", "#F39B7FFF" , "#00A087FF" ,"#8491B4FF","#CCCCCC","#4DBBD5FF","darkred","#B09C85FF" ,"#91D1C2FF"])'

sankeyNetwork(Links = data.frame(h3n2_jump3), Nodes = name,iterations = 0,
              colourScale = my_color,
              Source = "source", Target = "target", Value = "mean", NodeID = "name",fontFamily = NULL,
              sinksRight=FALSE, nodeWidth=10, nodePadding=5,height = 400, fontSize = 0) 

#sankey plot (epoch 2)
data_tmp <- h3n2_jump2[h3n2_jump2$period == "Period 2",] %>% group_by(startLocation) %>% summarise(n = sum(mean))
data_tmp1 <- h3n2_jump2[h3n2_jump2$period == "Period 2",] %>% group_by(endLocation) %>% summarise(n = sum(mean))
name_epoch2 <- data.frame(name = c(data_tmp$startLocation[order(data_tmp$n, decreasing = T)],
                                   data_tmp1$endLocation[order(data_tmp1$n, decreasing = T)]))
names1 <- data.frame(location = name_epoch2$name[1:12], source = 0:11)
names2 <- data.frame(location = name_epoch2$name[13:24], target = 12:23)
h3n2_jump4 <- left_join(h3n2_jump2[h3n2_jump2$period == "Period 2",], names1, by = c("startLocation" = "location"))
h3n2_jump4 <- left_join(h3n2_jump4, names2, by = c("endLocation" = "location"))
name_epoch2$name <- str_remove_all(name_epoch2$name," |/|-")

sankeyNetwork(Links = data.frame(h3n2_jump4), Nodes = name_epoch2,iterations = 0,
              colourScale = my_color,
              Source = "source", Target = "target", Value = "mean", NodeID = "name",fontFamily = NULL,
              sinksRight=FALSE, nodeWidth=10, nodePadding=5,height = 400, fontSize = 0) 

#sankey plot (epoch 3)
data_tmp <- h3n2_jump2[h3n2_jump2$period == "Period 3",] %>% group_by(startLocation) %>% summarise(n = sum(mean))
data_tmp1 <- h3n2_jump2[h3n2_jump2$period == "Period 3",] %>% group_by(endLocation) %>% summarise(n = sum(mean))
name_epoch3 <- data.frame(name = c(data_tmp$startLocation[order(data_tmp$n, decreasing = T)],
                                   data_tmp1$endLocation[order(data_tmp1$n, decreasing = T)]))
names1 <- data.frame(location = name_epoch3$name[1:12], source = 0:11)
names2 <- data.frame(location = name_epoch3$name[13:24], target = 12:23)
h3n2_jump5 <- left_join(h3n2_jump2[h3n2_jump2$period == "Period 3",], names1, by = c("startLocation" = "location"))
h3n2_jump5 <- left_join(h3n2_jump5, names2, by = c("endLocation" = "location"))
name_epoch3$name <- str_remove_all(name_epoch3$name," |/|-")

sankeyNetwork(Links = data.frame(h3n2_jump5), Nodes = name_epoch3,iterations = 0,
              colourScale = my_color,
              Source = "source", Target = "target", Value = "mean", NodeID = "name",fontFamily = NULL,
              sinksRight=FALSE, nodeWidth=10, nodePadding=5,height = 400, fontSize = 0) 

#==Net export plot==
h3n2_jump3 <- h3n2_jump1 %>%
  group_by(treeId, period, startLocation) %>%
  summarise(export_n = sum(n))
h3n2_jump3_1 <- h3n2_jump1 %>%
  group_by(treeId, period, endLocation) %>%
  summarise(import_n = sum(n))
h3n2_jump4 <- full_join(h3n2_jump3, h3n2_jump3_1, by = c("treeId" = "treeId", "period" = "period", "startLocation" = "endLocation"))
h3n2_jump4[is.na(h3n2_jump4)] <- 0
h3n2_jump5 <- h3n2_jump4 %>%
  mutate(net_export = export_n - import_n) %>%
  group_by(period, startLocation) %>%
  summarise(net_export_mean = mean(net_export),
            net_export_sd = sd(net_export),
            net_export_1q = quantile(net_export,0.25),
            net_export_3q = quantile(net_export,0.75))

h3n2_jump5 <- h3n2_jump5[order(h3n2_jump5$net_export_mean, decreasing = F),]
h3n2_jump5_1 <- h3n2_jump5[h3n2_jump5$period == "Period 1",]
factor(h3n2_jump5_1$startLocation, levels = h3n2_jump5_1$startLocation) -> h3n2_jump5_1$startLocation
ggplot(data = h3n2_jump5_1) +
  geom_hline(yintercept = 0, linetype = 2)+
  geom_bar(aes(x = startLocation, y = net_export_mean, fill = startLocation), stat = "identity", width = 0.5)+
  # geom_errorbar(aes(x = startLocation, ymin = net_export_mean-net_export_sd, ymax = net_export_mean+net_export_sd), width = 0.2)+
  geom_errorbar(aes(x = startLocation, ymin = net_export_1q, ymax = net_export_3q), width = 0.2, size = 0.2)+
  theme_bw()+
  theme(axis.text.y = element_blank(),
        plot.margin =  margin(0.1, 0.1, 0.1, 0.1, "cm"))+
  coord_flip()+
  scale_y_continuous(limits = c(min(h3n2_jump5$net_export_1q), max(h3n2_jump5$net_export_3q)), breaks = seq(-30,30,30))+
  scale_fill_manual(values = value)+
  labs(x = "", y = "Net exports") +
  guides(fill = F)-> fig4_1

h3n2_jump5_2 <- h3n2_jump5[h3n2_jump5$period == "Period 2",]
factor(h3n2_jump5_2$startLocation, levels = h3n2_jump5_2$startLocation) -> h3n2_jump5_2$startLocation
ggplot(data = h3n2_jump5_2) +
  geom_hline(yintercept = 0, linetype = 2)+
  geom_bar(aes(x = startLocation, y = net_export_mean, fill = startLocation), stat = "identity", width = 0.5)+
  geom_errorbar(aes(x = startLocation, ymin = net_export_1q, ymax = net_export_3q), width = 0.2, size = 0.2)+
  theme_bw()+
  theme(axis.text.y = element_blank(),
        plot.margin =  margin(0.1, 0.1, 0.1, 0.1, "cm"))+
  coord_flip()+
  scale_y_continuous(limits = c(min(h3n2_jump5$net_export_1q), max(h3n2_jump5$net_export_3q)), breaks = seq(-30,30,30))+
  scale_fill_manual(values = value)+
  labs(x = "", y = "Net exports") +
  guides(fill = F) -> fig5_1

h3n2_jump5_3 <- h3n2_jump5[h3n2_jump5$period == "Period 3",]
factor(h3n2_jump5_3$startLocation, levels = h3n2_jump5_3$startLocation) -> h3n2_jump5_3$startLocation
ggplot(data = h3n2_jump5_3) +
  geom_hline(yintercept = 0, linetype = 2)+
  geom_bar(aes(x = startLocation, y = net_export_mean, fill = startLocation), stat = "identity", width = 0.5)+
  geom_errorbar(aes(x = startLocation, ymin = net_export_1q, ymax = net_export_3q), width = 0.2, size = 0.2)+
  theme_bw()+
  theme(axis.text.y = element_blank(),
        plot.margin =  margin(0.1, 0.1, 0.1, 0.1, "cm"))+
  coord_flip()+
  scale_y_continuous(limits = c(min(h3n2_jump5$net_export_1q), max(h3n2_jump5$net_export_3q)), breaks = seq(-30,30,30))+
  scale_fill_manual(values = value)+
  labs(x = "", y = "Net exports") +
  guides(fill = F) -> fig6_1

#panel a (overall markov jumps)
h3n2_jump6 <- as.data.frame(h3n2_jump) %>%
  mutate(ISO_YEAR = isoyear(date),
         ISO_WEEK = isoweek(date),
         treeId = as.character(treeId)) %>%
  filter(!is.na(period)) %>%
  group_by(treeId, ISO_YEAR, ISO_WEEK) %>%
  summarise(n = n())
h3n2_jump6 <- h3n2_jump6 %>%
  left_join(date_match) %>%
  group_by(date_week) %>%
  summarise(mean = mean(n),
            sd = sd(n))

ggplot(data = h3n2_jump6) +
  annotate("rect", xmin = as.Date("2020-02-01"),xmax = as.Date("2021-07-31"),
           ymin = 0,ymax = 25,alpha = 0.2,fill = colors[5])+
  # geom_vline(xintercept = c(as.Date("2020-02-01"),as.Date("2021-08-01")), linetype = 2, size = 0.2)+
  geom_ribbon(aes(x = date_week, ymin = 0, ymax = mean), fill = colors1[2], alpha = 0.7)+
  geom_line(aes(x = date_week, y = mean, group = 1))+
  scale_y_continuous("Markov jumps", limits = c(-0.5,25), expand = c(0,0), breaks = seq(0,25,10))+
  scale_x_date(limits = c(as.Date("2018-02-01"), as.Date("2023-07-31")), expand = c(0.01,0))+
  theme_bw()+
  annotate("text", x = as.Date("2018-03-01"), y = 20, label = "H3N2", hjust = 0, size = 4)+
  theme(axis.title = element_blank(),
        plot.margin =  margin(0.1, 0.1, 0.1, 0.1, "cm"),
        axis.text.x = element_blank(),
        axis.text = element_text(size = 12),
        panel.grid.minor = element_blank())  -> m2

#==============
#==B/Victoria==
#==============
bv_jump <- read.delim("../genomic_part/post-analyses/jump_history/bv_even_Jumps.txt", sep = ",")
bv_jump$date <- decimal2Date(decimal_date(as.Date("2023-07-31")) - bv_jump$time)
bv_jump$period[bv_jump$date < as.Date("2020-02-01") & bv_jump$date >= as.Date("2018-02-01")] <- "Period 1"
bv_jump$period[bv_jump$date < as.Date("2021-08-01") & bv_jump$date >= as.Date("2020-02-01")] <- "Period 2"
bv_jump$period[bv_jump$date >= as.Date("2021-08-01")] <- "Period 3"
bv_jump1 <- bv_jump %>%
  filter(!is.na(period)) %>%
  group_by(treeId, period, startLocation, endLocation) %>%
  summarise(n = n())
bv_jump2 <- bv_jump1 %>%
  group_by(period, startLocation, endLocation) %>%
  summarise(mean = mean(n),
            sd = mean(n))
bv_jump2$startLocation[bv_jump2$startLocation == "SouthChina"] <- "China (S)"
bv_jump2$startLocation[bv_jump2$startLocation == "SoutheasternAsia"] <- "South-eastern Asia"
bv_jump2$startLocation[bv_jump2$startLocation == "SouthernAmerica"] <- "South Am"
bv_jump2$startLocation[bv_jump2$startLocation == "SouthernAsia"] <- "Southern Asia"
bv_jump2$startLocation[bv_jump2$startLocation == "WesternAsia"] <- "Western Asia"
bv_jump2$startLocation[bv_jump2$startLocation == "NorthChina"] <- "China (N)"
bv_jump2$startLocation[bv_jump2$startLocation == "NorthernAmerica"] <- "North Am"
bv_jump2$endLocation[bv_jump2$endLocation == "SouthChina"] <- "China (S)"
bv_jump2$endLocation[bv_jump2$endLocation == "SoutheasternAsia"] <- "South-eastern Asia"
bv_jump2$endLocation[bv_jump2$endLocation == "SouthernAmerica"] <- "South Am"
bv_jump2$endLocation[bv_jump2$endLocation == "SouthernAsia"] <- "Southern Asia"
bv_jump2$endLocation[bv_jump2$endLocation == "WesternAsia"] <- "Western Asia"
bv_jump2$endLocation[bv_jump2$endLocation == "NorthChina"] <- "China (N)"
bv_jump2$endLocation[bv_jump2$endLocation == "NorthernAmerica"] <- "North Am"

#sankey plot (epoch 1)
data_tmp <- bv_jump2[bv_jump2$period == "Period 1",] %>% group_by(startLocation) %>% summarise(n = sum(mean))
data_tmp1 <- bv_jump2[bv_jump2$period == "Period 1",] %>% group_by(endLocation) %>% summarise(n = sum(mean))
name <- data.frame(name = c(data_tmp$startLocation[order(data_tmp$n, decreasing = T)],
                            data_tmp1$endLocation[order(data_tmp1$n, decreasing = T)]))
names1 <- data.frame(location = name$name[1:12], source = 0:11)
names2 <- data.frame(location = name$name[13:24], target = 12:23)
bv_jump3 <- left_join(bv_jump2[bv_jump2$period == "Period 1",], names1, by = c("startLocation" = "location"))
bv_jump3 <- left_join(bv_jump3, names2, by = c("endLocation" = "location"))
name$name <- str_remove_all(name$name," |/|-")
my_color <- 'd3.scaleOrdinal() .domain(["China(N)","SoutheasternAsia","China(S)","JapanKorea", "SouthernAsia", "WesternAsia","NorthAm", "SouthAm",  "Europe",  "Russia", "Africa",  "Oceania"]) .range(["#7E6148FF","#3C5488FF","#FADDA9","#E64B35FF", "#F39B7FFF" , "#00A087FF" ,"#8491B4FF","#CCCCCC","#4DBBD5FF","darkred","#B09C85FF" ,"#91D1C2FF"])'

sankeyNetwork(Links = data.frame(bv_jump3), Nodes = name,iterations = 0,
              colourScale = my_color,
              Source = "source", Target = "target", Value = "mean", NodeID = "name",fontFamily = NULL,
              sinksRight=FALSE, nodeWidth=10, nodePadding=5,height = 400, fontSize = 0) 

#sankey plot (epoch 2)
data_tmp <- bv_jump2[bv_jump2$period == "Period 2",] %>% group_by(startLocation) %>% summarise(n = sum(mean))
data_tmp1 <- bv_jump2[bv_jump2$period == "Period 2",] %>% group_by(endLocation) %>% summarise(n = sum(mean))
name_epoch2 <- data.frame(name = c(data_tmp$startLocation[order(data_tmp$n, decreasing = T)],
                                   data_tmp1$endLocation[order(data_tmp1$n, decreasing = T)]))
names1 <- data.frame(location = name_epoch2$name[1:12], source = 0:11)
names2 <- data.frame(location = name_epoch2$name[13:24], target = 12:23)
bv_jump4 <- left_join(bv_jump2[bv_jump2$period == "Period 2",], names1, by = c("startLocation" = "location"))
bv_jump4 <- left_join(bv_jump4, names2, by = c("endLocation" = "location"))
name_epoch2$name <- str_remove_all(name_epoch2$name," |/|-")

sankeyNetwork(Links = data.frame(bv_jump4), Nodes = name_epoch2,iterations = 0,
              colourScale = my_color,
              Source = "source", Target = "target", Value = "mean", NodeID = "name",fontFamily = NULL,
              sinksRight=FALSE, nodeWidth=10, nodePadding=5,height = 400, fontSize = 0) 

#sankey plot (epoch 3)
data_tmp <- bv_jump2[bv_jump2$period == "Period 3",] %>% group_by(startLocation) %>% summarise(n = sum(mean))
data_tmp1 <- bv_jump2[bv_jump2$period == "Period 3",] %>% group_by(endLocation) %>% summarise(n = sum(mean))
name_epoch3 <- data.frame(name = c(data_tmp$startLocation[order(data_tmp$n, decreasing = T)],
                                   data_tmp1$endLocation[order(data_tmp1$n, decreasing = T)]))
names1 <- data.frame(location = name_epoch3$name[1:12], source = 0:11)
names2 <- data.frame(location = name_epoch3$name[13:24], target = 12:23)
bv_jump5 <- left_join(bv_jump2[bv_jump2$period == "Period 3",], names1, by = c("startLocation" = "location"))
bv_jump5 <- left_join(bv_jump5, names2, by = c("endLocation" = "location"))
name_epoch3$name <- str_remove_all(name_epoch3$name," |/|-")

sankeyNetwork(Links = data.frame(bv_jump5), Nodes = name_epoch3,iterations = 0,
              colourScale = my_color,
              Source = "source", Target = "target", Value = "mean", NodeID = "name",fontFamily = NULL,
              sinksRight=FALSE, nodeWidth=10, nodePadding=5,height = 400, fontSize = 0) 

#==Net export plot==
bv_jump3 <- bv_jump1 %>%
  group_by(treeId, period, startLocation) %>%
  summarise(export_n = sum(n))
bv_jump3_1 <- bv_jump1 %>%
  group_by(treeId, period, endLocation) %>%
  summarise(import_n = sum(n))
bv_jump4 <- full_join(bv_jump3, bv_jump3_1, by = c("treeId" = "treeId", "period" = "period", "startLocation" = "endLocation"))
bv_jump4[is.na(bv_jump4)] <- 0
bv_jump5 <- bv_jump4 %>%
  mutate(net_export = export_n - import_n) %>%
  group_by(period, startLocation) %>%
  summarise(net_export_mean = mean(net_export),
            net_export_sd = sd(net_export),
            net_export_1q = quantile(net_export,0.25),
            net_export_3q = quantile(net_export,0.75))

bv_jump5 <- bv_jump5[order(bv_jump5$net_export_mean, decreasing = F),]
bv_jump5_1 <- bv_jump5[bv_jump5$period == "Period 1",]
factor(bv_jump5_1$startLocation, levels = bv_jump5_1$startLocation) -> bv_jump5_1$startLocation
ggplot(data = bv_jump5_1) +
  geom_hline(yintercept = 0, linetype = 2)+
  geom_bar(aes(x = startLocation, y = net_export_mean, fill = startLocation), stat = "identity", width = 0.5)+
  # geom_errorbar(aes(x = startLocation, ymin = net_export_mean-net_export_sd, ymax = net_export_mean+net_export_sd), width = 0.2)+
  geom_errorbar(aes(x = startLocation, ymin = net_export_1q, ymax = net_export_3q), width = 0.2, size = 0.2)+
  theme_bw()+
  theme(axis.text.y = element_blank(),
        plot.margin =  margin(0.1, 0.1, 0.1, 0.1, "cm"))+
  coord_flip()+
  scale_y_continuous(limits = c(min(bv_jump5$net_export_1q), max(bv_jump5$net_export_3q)), breaks = seq(-30,30,30))+
  scale_fill_manual(values = value)+
  labs(x = "", y = "Net exports") +
  guides(fill = F)-> fig7_1

bv_jump5_2 <- bv_jump5[bv_jump5$period == "Period 2",]
factor(bv_jump5_2$startLocation, levels = bv_jump5_2$startLocation) -> bv_jump5_2$startLocation
ggplot(data = bv_jump5_2) +
  geom_hline(yintercept = 0, linetype = 2)+
  geom_bar(aes(x = startLocation, y = net_export_mean, fill = startLocation), stat = "identity", width = 0.5)+
  geom_errorbar(aes(x = startLocation, ymin = net_export_1q, ymax = net_export_3q), width = 0.2, size = 0.2)+
  theme_bw()+
  theme(axis.text.y = element_blank(),
        plot.margin =  margin(0.1, 0.1, 0.1, 0.1, "cm"))+
  coord_flip()+
  scale_y_continuous(limits = c(min(bv_jump5$net_export_1q), max(bv_jump5$net_export_3q)), breaks = seq(-30,30,30))+
  scale_fill_manual(values = value)+
  labs(x = "", y = "Net exports") +
  guides(fill = F) -> fig8_1

bv_jump5_3 <- bv_jump5[bv_jump5$period == "Period 3",]
factor(bv_jump5_3$startLocation, levels = bv_jump5_3$startLocation) -> bv_jump5_3$startLocation
ggplot(data = bv_jump5_3) +
  geom_hline(yintercept = 0, linetype = 2)+
  geom_bar(aes(x = startLocation, y = net_export_mean, fill = startLocation), stat = "identity", width = 0.5)+
  geom_errorbar(aes(x = startLocation, ymin = net_export_1q, ymax = net_export_3q), width = 0.2, size = 0.2)+
  theme_bw()+
  theme(axis.text.y = element_blank(),
        plot.margin =  margin(0.1, 0.1, 0.1, 0.1, "cm"))+
  coord_flip()+
  scale_y_continuous(limits = c(min(bv_jump5$net_export_1q), max(bv_jump5$net_export_3q)), breaks = seq(-30,30,30))+
  scale_fill_manual(values = value)+
  labs(x = "", y = "Net exports") +
  guides(fill = F) -> fig9_1

#panel a (overall markov jumps)
bv_jump6 <- as.data.frame(bv_jump) %>%
  mutate(ISO_YEAR = isoyear(date),
         ISO_WEEK = isoweek(date),
         treeId = as.character(treeId)) %>%
  filter(!is.na(period)) %>%
  group_by(treeId, ISO_YEAR, ISO_WEEK) %>%
  summarise(n = n())
bv_jump6 <- bv_jump6 %>%
  left_join(date_match) %>%
  group_by(date_week) %>%
  summarise(mean = mean(n),
            sd = sd(n))

ggplot(data = bv_jump6) +
  annotate("rect", xmin = as.Date("2020-02-01"),xmax = as.Date("2021-07-31"),
           ymin = 0,ymax = 25,alpha = 0.2,fill = colors[5])+
  # geom_vline(xintercept = c(as.Date("2020-02-01"),as.Date("2021-08-01")), linetype = 2, size = 0.2)+
  geom_ribbon(aes(x = date_week, ymin = 0, ymax = mean), fill = colors1[3], alpha = 0.7)+
  geom_line(aes(x = date_week, y = mean, group = 1))+
  scale_y_continuous("Markov jumps", limits = c(-0.5,25), expand = c(0,0), breaks = seq(0,25,10))+
  scale_x_date(limits = c(as.Date("2018-02-01"), as.Date("2023-07-31")), expand = c(0.01,0))+
  theme_bw()+
  annotate("text", x = as.Date("2018-03-01"), y = 20, label = "B/Victoria", hjust = 0, size = 4)+
  theme(axis.title = element_blank(),
        plot.margin =  margin(0.1, 0.1, 0.1, 0.1, "cm"),
        axis.text.x = element_blank(),
        axis.text = element_text(size = 12),
        panel.grid.minor = element_blank())  -> m3

#==============
#==B/Yamagata==
#==============
by_jump <- read.delim("../genomic_part/post-analyses/jump_history/by_even_Jumps.txt", sep = ",")
by_jump$date <- decimal2Date(decimal_date(as.Date("2020-03-24")) - by_jump$time)
by_jump$period[by_jump$date < as.Date("2020-02-01") & by_jump$date >= as.Date("2018-02-01")] <- "Period 1"
by_jump$period[by_jump$date < as.Date("2021-08-01") & by_jump$date >= as.Date("2020-02-01")] <- "Period 2"
by_jump$period[by_jump$date >= as.Date("2021-08-01")] <- "Period 3"
by_jump1 <- by_jump %>%
  filter(!is.na(period)) %>%
  group_by(treeId, period, startLocation, endLocation) %>%
  summarise(n = n())
by_jump2 <- by_jump1 %>%
  group_by(period, startLocation, endLocation) %>%
  summarise(mean = mean(n),
            sd = mean(n))
by_jump2$startLocation[by_jump2$startLocation == "SouthChina"] <- "China (S)"
by_jump2$startLocation[by_jump2$startLocation == "SoutheasternAsia"] <- "South-eastern Asia"
by_jump2$startLocation[by_jump2$startLocation == "SouthernAmerica"] <- "South Am"
by_jump2$startLocation[by_jump2$startLocation == "SouthernAsia"] <- "Southern Asia"
by_jump2$startLocation[by_jump2$startLocation == "WesternAsia"] <- "Western Asia"
by_jump2$startLocation[by_jump2$startLocation == "NorthChina"] <- "China (N)"
by_jump2$startLocation[by_jump2$startLocation == "NorthernAmerica"] <- "North Am"
by_jump2$endLocation[by_jump2$endLocation == "SouthChina"] <- "China (S)"
by_jump2$endLocation[by_jump2$endLocation == "SoutheasternAsia"] <- "South-eastern Asia"
by_jump2$endLocation[by_jump2$endLocation == "SouthernAmerica"] <- "South Am"
by_jump2$endLocation[by_jump2$endLocation == "SouthernAsia"] <- "Southern Asia"
by_jump2$endLocation[by_jump2$endLocation == "WesternAsia"] <- "Western Asia"
by_jump2$endLocation[by_jump2$endLocation == "NorthChina"] <- "China (N)"
by_jump2$endLocation[by_jump2$endLocation == "NorthernAmerica"] <- "North Am"

#sankey plot (epoch 1)
data_tmp <- by_jump2[by_jump2$period == "Period 1",] %>% group_by(startLocation) %>% summarise(n = sum(mean))
data_tmp1 <- by_jump2[by_jump2$period == "Period 1",] %>% group_by(endLocation) %>% summarise(n = sum(mean))
name <- data.frame(name = c(data_tmp$startLocation[order(data_tmp$n, decreasing = T)],
                            data_tmp1$endLocation[order(data_tmp1$n, decreasing = T)]))
names1 <- data.frame(location = name$name[1:12], source = 0:11)
names2 <- data.frame(location = name$name[13:24], target = 12:23)
by_jump3 <- left_join(by_jump2[by_jump2$period == "Period 1",], names1, by = c("startLocation" = "location"))
by_jump3 <- left_join(by_jump3, names2, by = c("endLocation" = "location"))
name$name <- str_remove_all(name$name," |/|-")
my_color <- 'd3.scaleOrdinal() .domain(["China(N)","SoutheasternAsia","China(S)","JapanKorea", "SouthernAsia", "WesternAsia","NorthAm", "SouthAm",  "Europe",  "Russia", "Africa",  "Oceania"]) .range(["#7E6148FF","#3C5488FF","#FADDA9","#E64B35FF", "#F39B7FFF" , "#00A087FF" ,"#8491B4FF","#CCCCCC","#4DBBD5FF","darkred","#B09C85FF" ,"#91D1C2FF"])'

sankeyNetwork(Links = data.frame(by_jump3), Nodes = name,iterations = 0,
              colourScale = my_color,
              Source = "source", Target = "target", Value = "mean", NodeID = "name",fontFamily = NULL,
              sinksRight=FALSE, nodeWidth=10, nodePadding=5, height = 400,fontSize = 0) 

#==Net export plot==
by_jump3 <- by_jump1 %>%
  group_by(treeId, period, startLocation) %>%
  summarise(export_n = sum(n))
by_jump3_1 <- by_jump1 %>%
  group_by(treeId, period, endLocation) %>%
  summarise(import_n = sum(n))
by_jump4 <- full_join(by_jump3, by_jump3_1, by = c("treeId" = "treeId", "period" = "period", "startLocation" = "endLocation"))
by_jump4[is.na(by_jump4)] <- 0
by_jump5 <- by_jump4 %>%
  mutate(net_export = export_n - import_n) %>%
  group_by(period, startLocation) %>%
  summarise(net_export_mean = mean(net_export),
            net_export_sd = sd(net_export),
            net_export_1q = quantile(net_export,0.25),
            net_export_3q = quantile(net_export,0.75))

by_jump5 <- by_jump5[order(by_jump5$net_export_mean, decreasing = F),]
by_jump5_1 <- by_jump5[by_jump5$period == "Period 1",]
factor(by_jump5_1$startLocation, levels = by_jump5_1$startLocation) -> by_jump5_1$startLocation
ggplot(data = by_jump5_1) +
  geom_hline(yintercept = 0, linetype = 2)+
  geom_bar(aes(x = startLocation, y = net_export_mean, fill = startLocation), stat = "identity", width = 0.5)+
  # geom_errorbar(aes(x = startLocation, ymin = net_export_mean-net_export_sd, ymax = net_export_mean+net_export_sd), width = 0.2)+
  geom_errorbar(aes(x = startLocation, ymin = net_export_1q, ymax = net_export_3q), width = 0.2, size = 0.2)+
  theme_bw()+
  theme(axis.text.y = element_blank(),
        plot.margin =  margin(0.1, 0.1, 0.1, 0.1, "cm"))+
  coord_flip()+
  scale_y_continuous(limits = c(min(by_jump5$net_export_1q), max(by_jump5$net_export_3q)), breaks = seq(-30,30,30))+
  scale_fill_manual(values = value)+
  labs(x = "", y = "Net exports") +
  guides(fill = F)-> fig10_1

#panel a (overall markov jumps)
by_jump6 <- as.data.frame(by_jump) %>%
  mutate(ISO_YEAR = isoyear(date),
         ISO_WEEK = isoweek(date),
         treeId = as.character(treeId)) %>%
  filter(!is.na(period)) %>%
  group_by(treeId, ISO_YEAR, ISO_WEEK) %>%
  summarise(n = n())
by_jump6 <- by_jump6 %>%
  left_join(date_match) %>%
  group_by(date_week) %>%
  summarise(mean = mean(n),
            sd = sd(n))

ggplot(data = by_jump6) +
  annotate("rect", xmin = as.Date("2020-02-01"),xmax = as.Date("2021-07-31"),
           ymin = 0,ymax = 25,alpha = 0.2,fill = colors[5])+
  # geom_vline(xintercept = c(as.Date("2020-02-01"),as.Date("2021-08-01")), linetype = 2, size = 0.2)+
  geom_ribbon(aes(x = date_week, ymin = 0, ymax = mean), fill = colors1[9], alpha = 0.7)+
  geom_line(aes(x = date_week, y = mean, group = 1))+
  scale_y_continuous("Markov jumps", limits = c(-0.5,25), expand = c(0,0), breaks = seq(0,25,10))+
  scale_x_date(limits = c(as.Date("2018-02-01"), as.Date("2023-07-31")), expand = c(0.01,0))+
  theme_bw()+
  annotate("text", x = as.Date("2018-03-01"), y = 20, label = "B/Yamagata", hjust = 0, size = 4)+
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 12),
        plot.margin =  margin(0.1, 0.1, 0.1, 0.1, "cm"),
        panel.grid.minor = element_blank())  -> m4

#==output==
svg("output/Fig4.svg",width = 10,height = 14)
(fig_white|fig1_1|fig_white|fig2_1|fig_white|fig3_1)/
  (fig_white|fig4_1|fig_white|fig5_1|fig_white|fig6_1)/
  (fig_white|fig7_1|fig_white|fig8_1|fig_white|fig9_1)/
  (fig_white|fig10_1|fig_white|fig_white|fig_white|fig_white) -> part2
(m1/m2/m3/m4)+plot_annotation(subtitle = "a. Counts of Markov Jumps")&theme(plot.subtitle = element_text(size = 12, face = "bold")) -> part1
# viewport(x = 0, y = 0.78, width = 1, height = 0.22, just = c("left", "bottom")) -> vp1
# viewport(x = 0, y = 0, width = 1, height = 0.8, just = c("left", "bottom")) -> vp2
# viewport(x = 0.35, y = 0.01, width = 0.6, height = 0.2, just = c("left", "bottom")) -> vp3
viewport(x = 0, y = 0.72, width = 1, height = 0.28, just = c("left", "bottom")) -> vp1
viewport(x = 0, y = 0, width = 1, height = 0.74, just = c("left", "bottom")) -> vp2
viewport(x = 0.35, y = 0.008, width = 0.6, height = 0.2, just = c("left", "bottom")) -> vp3
print(part2,vp = vp2)
print(part1,vp = vp1)
print(fig_legend, vp = vp3)
dev.off()
