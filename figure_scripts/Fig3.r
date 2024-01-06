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

#==h1n1pdm09==
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
h1n1_jump2$startLocation[h1n1_jump2$startLocation == "SouthChina"] <- "South China"
h1n1_jump2$startLocation[h1n1_jump2$startLocation == "SoutheasternAsia"] <- "Southeast Asia"
h1n1_jump2$startLocation[h1n1_jump2$startLocation == "SouthernAmerica"] <- "South America"
h1n1_jump2$startLocation[h1n1_jump2$startLocation == "SouthernAsia"] <- "South Asia"
h1n1_jump2$startLocation[h1n1_jump2$startLocation == "WesternAsia"] <- "West Asia"
h1n1_jump2$startLocation[h1n1_jump2$startLocation == "NorthChina"] <- "North China"
h1n1_jump2$startLocation[h1n1_jump2$startLocation == "NorthernAmerica"] <- "North America"
h1n1_jump2$endLocation[h1n1_jump2$endLocation == "SouthChina"] <- "South China"
h1n1_jump2$endLocation[h1n1_jump2$endLocation == "SoutheasternAsia"] <- "Southeast Asia"
h1n1_jump2$endLocation[h1n1_jump2$endLocation == "SouthernAmerica"] <- "South America"
h1n1_jump2$endLocation[h1n1_jump2$endLocation == "SouthernAsia"] <- "South Asia"
h1n1_jump2$endLocation[h1n1_jump2$endLocation == "WesternAsia"] <- "West Asia"
h1n1_jump2$endLocation[h1n1_jump2$endLocation == "NorthChina"] <- "North China"
h1n1_jump2$endLocation[h1n1_jump2$endLocation == "NorthernAmerica"] <- "North America"
levels <- c("Japan/Korea","North China","South China", "South Asia", "West Asia","Southeast Asia","Oceania",
            "North America", "South America",  "Europe",  "Russia", "Africa")

##
factor(h1n1_jump2$startLocation, levels = levels) -> h1n1_jump2$startLocation
factor(h1n1_jump2$endLocation, levels = levels) -> h1n1_jump2$endLocation
cut(h1n1_jump2$mean, breaks = c(0,2.5,5,10,20,30,45),right = T,
    labels = c("(0, 2.5]","(2.5, 5]","(5, 10]","(10, 20]","(20, 30]","(30, 45]")) -> h1n1_jump2$mean_1

ggplot()+
  geom_tile(data = h1n1_jump2[h1n1_jump2$period == "Period 1",],
            aes(x = endLocation, y = startLocation, fill = mean_1))+
  # scale_fill_gradientn("Count of transition event",colors = rev(c("#E64B35FF","#0061A3", "#CDDDE5","grey93")),
  #                      na.value="white",limits = c(1, 42))+
  scale_fill_manual("Count of transition event",values = c("grey90","#CDDDE5", "#7BABCA", "#2879B0", "#2E5C8C",  "#E64B35"),
                    na.value="white")+
  scale_x_discrete("To")+
  scale_y_discrete("H1N1pdm09\n\nFrom")+
  theme_bw()+
  # guides(fill = "none")+
  guides(fill = guide_colorbar(frame.colour = "black",
                               ticks.colour = "black"))+
  labs(title = "Period 1 (Feb 2018 to Jan 2020)")+
  theme(
    # axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
    legend.key.width = unit(1,"cm"),
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none")+
  labs(tag = "b")-> p1

ggplot()+
  geom_tile(data = h1n1_jump2[h1n1_jump2$period == "Period 2",],
            aes(x = endLocation, y = startLocation, fill = mean_1))+
  # scale_fill_gradientn("Count of transition event",colors = rev(c("#E64B35FF","#0061A3", "#CDDDE5","grey93")),
  #                      na.value="white",limits = c(1, 42))+
  scale_fill_manual("Count of transition event",values = c("grey90","#CDDDE5", "#7BABCA", "#2879B0", "#2E5C8C",  "#E64B35"),
                    na.value="white")+
  scale_x_discrete("To")+
  scale_y_discrete("From")+
  theme_bw()+
  # guides(fill = "none")+
  guides(fill = guide_colorbar(frame.colour = "black",
                               ticks.colour = "black"))+
  labs(title = "Period 2 (Feb 2020 to Jul 2021)")+
  theme(
    # axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
    legend.key.width = unit(1,"cm"),
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none") -> p2

ggplot()+
  geom_tile(data = h1n1_jump2[h1n1_jump2$period == "Period 3",],
            aes(x = endLocation, y = startLocation, fill = mean_1))+
  # scale_fill_gradientn("Count of transition event",colors = rev(c("#E64B35FF","#0061A3", "#CDDDE5","grey93")),
  #                      na.value="white",limits = c(1, 42))+
  scale_fill_manual("Count of transition event",values = c("grey90","#CDDDE5", "#7BABCA", "#2879B0", "#E64B35"),
                    na.value="white")+
  scale_x_discrete("To")+
  scale_y_discrete("From")+
  theme_bw()+
  # guides(fill = "none")+
  guides(fill = guide_colorbar(frame.colour = "black",
                               ticks.colour = "black"))+
  labs(title = "Period 3 (Aug 2021 to Jul 2023)")+
  theme(
    # axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
    legend.key.width = unit(1,"cm"),
    plot.title = element_text(hjust = 0.5),
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none") -> p3
##
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

#panel a
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

#==h3n2==
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
h3n2_jump2$startLocation[h3n2_jump2$startLocation == "SouthChina"] <- "South China"
h3n2_jump2$startLocation[h3n2_jump2$startLocation == "SoutheasternAsia"] <- "Southeast Asia"
h3n2_jump2$startLocation[h3n2_jump2$startLocation == "SouthernAmerica"] <- "South America"
h3n2_jump2$startLocation[h3n2_jump2$startLocation == "SouthernAsia"] <- "South Asia"
h3n2_jump2$startLocation[h3n2_jump2$startLocation == "WesternAsia"] <- "West Asia"
h3n2_jump2$startLocation[h3n2_jump2$startLocation == "NorthChina"] <- "North China"
h3n2_jump2$startLocation[h3n2_jump2$startLocation == "NorthernAmerica"] <- "North America"
h3n2_jump2$endLocation[h3n2_jump2$endLocation == "SouthChina"] <- "South China"
h3n2_jump2$endLocation[h3n2_jump2$endLocation == "SoutheasternAsia"] <- "Southeast Asia"
h3n2_jump2$endLocation[h3n2_jump2$endLocation == "SouthernAmerica"] <- "South America"
h3n2_jump2$endLocation[h3n2_jump2$endLocation == "SouthernAsia"] <- "South Asia"
h3n2_jump2$endLocation[h3n2_jump2$endLocation == "WesternAsia"] <- "West Asia"
h3n2_jump2$endLocation[h3n2_jump2$endLocation == "NorthChina"] <- "North China"
h3n2_jump2$endLocation[h3n2_jump2$endLocation == "NorthernAmerica"] <- "North America"

##
factor(h3n2_jump2$startLocation, levels = levels) -> h3n2_jump2$startLocation
factor(h3n2_jump2$endLocation, levels = levels) -> h3n2_jump2$endLocation
cut(h3n2_jump2$mean, breaks = c(0,2.5,5,10,20,30,45),right = T,
    labels = c("(0, 2.5]","(2.5, 5]","(5, 10]","(10, 20]","(20, 30]","(30, 45]")) -> h3n2_jump2$mean_1

ggplot()+
  geom_tile(data = h3n2_jump2[h3n2_jump2$period == "Period 1",],
            aes(x = endLocation, y = startLocation, fill = mean_1))+
  # scale_fill_gradientn("Count of transition event",colors = rev(c("#E64B35FF","#0061A3", "#CDDDE5","grey93")),
  #                      na.value="white",limits = c(1, 42))+
  scale_fill_manual("Count of transition event",values = c("grey90","#CDDDE5", "#7BABCA", "#2879B0", "#2E5C8C",  "#E64B35"),
                    na.value="white")+
  scale_x_discrete("To")+
  scale_y_discrete("H3N2\n\nFrom")+
  theme_bw()+
  # guides(fill = "none")+
  guides(fill = guide_colorbar(frame.colour = "black",
                               ticks.colour = "black"))+
  # labs(title = "Epoch 1 (Feb 2018 to Jan 2020)")+
  theme(
    # axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
    legend.key.width = unit(1,"cm"),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none")+
  labs(tag = "c") -> p4

ggplot()+
  geom_tile(data = h3n2_jump2[h3n2_jump2$period == "Period 2",],
            aes(x = endLocation, y = startLocation, fill = mean_1))+
  # scale_fill_gradientn("Count of transition event",colors = rev(c("#E64B35FF","#0061A3", "#CDDDE5","grey93")),
  #                      na.value="white",limits = c(1, 42))+
  scale_fill_manual("Count of transition event",values = c("grey90","#CDDDE5", "#7BABCA", "#2879B0", "#2E5C8C",  "#E64B35"),
                    na.value="white")+
  scale_x_discrete("To")+
  scale_y_discrete("From")+
  theme_bw()+
  # guides(fill = "none")+
  guides(fill = guide_colorbar(frame.colour = "black",
                               ticks.colour = "black"))+
  # labs(title = "Epoch 2 (Feb 2020 to Jul 2021)")+
  theme(
    # axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
    legend.key.width = unit(1,"cm"),
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none") -> p5

ggplot()+
  geom_tile(data = h3n2_jump2[h3n2_jump2$period == "Period 3",],
            aes(x = endLocation, y = startLocation, fill = mean_1))+
  # scale_fill_gradientn("Count of transition event",colors = rev(c("#E64B35FF","#0061A3", "#CDDDE5","grey93")),
  #                      na.value="white",limits = c(1, 42))+
  scale_fill_manual("Count of transition event",values = c("grey90","#CDDDE5", "#7BABCA", "#2879B0", "#2E5C8C",  "#E64B35"),
                    na.value="white")+
  scale_x_discrete("To")+
  scale_y_discrete("From")+
  theme_bw()+
  # guides(fill = "none")+
  guides(fill = guide_colorbar(frame.colour = "black",
                               ticks.colour = "black"))+
  # labs(title = "Epoch 3 (Aug 2021 to Jul 2023)")+
  theme(
    # axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
    legend.key.width = unit(1,"cm"),
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none") -> p6
##

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

#==B/Victoria==
bv_jump <- read.delim("../genomic_part/post-analyses/jump_history/bv_even_Jumps.txt", sep = ",")
bv_jump$date <- decimal2Date(decimal_date(as.Date("2023-07-31")) - bv_jump$time)
bv_jump$period[bv_jump$date < as.Date("2020-02-01") & bv_jump$date >= as.Date("2018-02-01")] <- "Period 1"
bv_jump$period[bv_jump$date < as.Date("2021-08-01") & bv_jump$date >= as.Date("2020-02-01")] <- "Period 2"
bv_jump$period[bv_jump$date >= as.Date("2021-08-01")] <- "Period 3"
# nrow(bv_jump[bv_jump$date < as.Date("2018-02-01"),])
bv_jump1 <- bv_jump %>%
  filter(!is.na(period)) %>%
  group_by(treeId, period, startLocation, endLocation) %>%
  summarise(n = n())
bv_jump2 <- bv_jump1 %>%
  group_by(period, startLocation, endLocation) %>%
  summarise(mean = mean(n),
            sd = mean(n))
bv_jump2$startLocation[bv_jump2$startLocation == "SouthChina"] <- "South China"
bv_jump2$startLocation[bv_jump2$startLocation == "SoutheasternAsia"] <- "Southeast Asia"
bv_jump2$startLocation[bv_jump2$startLocation == "SouthernAmerica"] <- "South America"
bv_jump2$startLocation[bv_jump2$startLocation == "SouthernAsia"] <- "South Asia"
bv_jump2$startLocation[bv_jump2$startLocation == "WesternAsia"] <- "West Asia"
bv_jump2$startLocation[bv_jump2$startLocation == "NorthChina"] <- "North China"
bv_jump2$startLocation[bv_jump2$startLocation == "NorthernAmerica"] <- "North America"
bv_jump2$endLocation[bv_jump2$endLocation == "SouthChina"] <- "South China"
bv_jump2$endLocation[bv_jump2$endLocation == "SoutheasternAsia"] <- "Southeast Asia"
bv_jump2$endLocation[bv_jump2$endLocation == "SouthernAmerica"] <- "South America"
bv_jump2$endLocation[bv_jump2$endLocation == "SouthernAsia"] <- "South Asia"
bv_jump2$endLocation[bv_jump2$endLocation == "WesternAsia"] <- "West Asia"
bv_jump2$endLocation[bv_jump2$endLocation == "NorthChina"] <- "North China"
bv_jump2$endLocation[bv_jump2$endLocation == "NorthernAmerica"] <- "North America"

##
factor(bv_jump2$startLocation, levels = levels) -> bv_jump2$startLocation
factor(bv_jump2$endLocation, levels = levels) -> bv_jump2$endLocation
cut(bv_jump2$mean, breaks = c(0,2.5,5,10,20,30,45),right = T,
    labels = c("(0, 2.5]","(2.5, 5]","(5, 10]","(10, 20]","(20, 30]","(30, 45]")) -> bv_jump2$mean_1

ggplot()+
  geom_tile(data = bv_jump2[bv_jump2$period == "Period 1",],
            aes(x = endLocation, y = startLocation, fill = mean_1))+
  # scale_fill_gradientn("Count of transition event",colors = rev(c("#E64B35FF","#0061A3", "#CDDDE5","grey93")),
  #                      na.value="white",limits = c(1, 42))+
  scale_fill_manual("Count of transition event",values = c("grey90","#CDDDE5", "#7BABCA", "#2879B0", "#2E5C8C",  "#E64B35"),
                    na.value="white")+
  scale_x_discrete("To")+
  scale_y_discrete("B/Victoria\n\nFrom")+
  theme_bw()+
  # guides(fill = "none")+
  guides(fill = guide_colorbar(frame.colour = "black",
                               ticks.colour = "black"))+
  # labs(title = "Epoch 1 (Feb 2018 to Jan 2020)")+
  theme(
    # axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
    legend.key.width = unit(1,"cm"),
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "none")+
  labs(tag = "d") -> p7

ggplot()+
  geom_tile(data = bv_jump2[bv_jump2$period == "Period 2",],
            aes(x = endLocation, y = startLocation, fill = mean_1))+
  # scale_fill_gradientn("Count of transition event",colors = rev(c("#E64B35FF","#0061A3", "#CDDDE5","grey93")),
  #                      na.value="white",limits = c(1, 42))+
  scale_fill_manual("Count of transition event",values = c("grey90","#CDDDE5", "#7BABCA", "#2879B0", "#2E5C8C",  "#E64B35"),
                    na.value="white")+
  scale_x_discrete("To")+
  scale_y_discrete("From")+
  theme_bw()+
  # guides(fill = "none")+
  guides(fill = guide_colorbar(frame.colour = "black",
                               ticks.colour = "black"))+
  # labs(title = "Epoch 2 (Feb 2020 to Jul 2021)")+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0),
        legend.key.width = unit(1,"cm"),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none") -> p8

ggplot()+
  geom_tile(data = bv_jump2[bv_jump2$period == "Period 3",],
            aes(x = endLocation, y = startLocation, fill = mean_1))+
  # scale_fill_gradientn("Count of transition event",colors = rev(c("#E64B35FF","#0061A3", "#CDDDE5","grey93")),
  #                      na.value="white",limits = c(1, 42))+
  scale_fill_manual("Count of transition event",values = c("grey90","#CDDDE5", "#7BABCA", "#2879B0", "#2E5C8C",  "#E64B35"),
                    na.value="white")+
  scale_x_discrete("To")+
  scale_y_discrete("From")+
  theme_bw()+
  # guides(fill = "none")+
  guides(fill = guide_colorbar(frame.colour = "black",
                               ticks.colour = "black"))+
  # labs(title = "Epoch 3 (Aug 2021 to Jul 2023)")+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0),
        legend.key.width = unit(1,"cm"),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none") -> p9

##
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

#==B/Yamagata==
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
by_jump2$startLocation[by_jump2$startLocation == "SouthChina"] <- "South China"
by_jump2$startLocation[by_jump2$startLocation == "SoutheasternAsia"] <- "Southeast Asia"
by_jump2$startLocation[by_jump2$startLocation == "SouthernAmerica"] <- "South America"
by_jump2$startLocation[by_jump2$startLocation == "SouthernAsia"] <- "South Asia"
by_jump2$startLocation[by_jump2$startLocation == "WesternAsia"] <- "West Asia"
by_jump2$startLocation[by_jump2$startLocation == "NorthChina"] <- "North China"
by_jump2$startLocation[by_jump2$startLocation == "NorthernAmerica"] <- "North America"
by_jump2$endLocation[by_jump2$endLocation == "SouthChina"] <- "South China"
by_jump2$endLocation[by_jump2$endLocation == "SoutheasternAsia"] <- "Southeast Asia"
by_jump2$endLocation[by_jump2$endLocation == "SouthernAmerica"] <- "South America"
by_jump2$endLocation[by_jump2$endLocation == "SouthernAsia"] <- "South Asia"
by_jump2$endLocation[by_jump2$endLocation == "WesternAsia"] <- "West Asia"
by_jump2$endLocation[by_jump2$endLocation == "NorthChina"] <- "North China"
by_jump2$endLocation[by_jump2$endLocation == "NorthernAmerica"] <- "North America"

##
factor(by_jump2$startLocation, levels = levels) -> by_jump2$startLocation
factor(by_jump2$endLocation, levels = levels) -> by_jump2$endLocation
cut(by_jump2$mean, breaks = c(0,2.5,5,10,20,30,45),right = T,
    labels = c("(0, 2.5]","(2.5, 5]","(5, 10]","(10, 20]","(20, 30]","(30, 45]")) -> by_jump2$mean_1
ggplot()+
  geom_tile(data = by_jump2[by_jump2$period == "Period 1",],
            aes(x = endLocation, y = startLocation, fill = mean_1))+
  # scale_fill_gradientn("Count of transition event",colors = rev(c("#E64B35FF","#0061A3", "#CDDDE5","grey93")),
  #                      na.value="white",limits = c(1, 42))+
  scale_fill_manual("Count of transition event",values = c("grey90","#CDDDE5", "#7BABCA", "#2879B0", "#2E5C8C",  "#E64B35"),
                    na.value="white")+
  scale_x_discrete("To")+
  scale_y_discrete("B/Yamagata\n\nFrom")+
  theme_bw()+
  # guides(fill = "none")+
  guides(fill = guide_colorbar(frame.colour = "black",
                               ticks.colour = "black",
                               title.position = "top"))+
  # labs(title = "Epoch 1 (Feb 2018 to Jan 2020)")+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0),
        legend.key.width = unit(1,"cm"),
        panel.grid = element_blank(),
        legend.direction = "horizontal",
        legend.position = "none")+
  labs(tag = "e")-> p10

ggplot()+
  geom_tile(data = h1n1_jump2[h1n1_jump2$period == "Period 1",],
            aes(x = endLocation, y = startLocation, fill = mean_1))+
  # scale_fill_gradientn("Count of transition event",colors = rev(c("#E64B35FF","#0061A3", "#CDDDE5","grey93")),
  #                      na.value="white",limits = c(1, 42))+
  scale_fill_manual("Counts of transition events",values = c("grey90","#CDDDE5", "#7BABCA", "#2879B0", "#2E5C8C",  "#E64B35"),
                    na.value="white")+
  scale_x_discrete("To")+
  scale_y_discrete("B/Yamagata\n\nFrom")+
  # theme_bw()+
  # guides(fill = "none")+
  guides(fill = guide_legend(
    # frame.colour = "black",
    # ticks.colour = "black",
    title.position = "top"))+
  # labs(title = "Epoch 1 (Feb 2018 to Jan 2020)")+
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    plot.background = element_rect(fill = "white",color = "transparent"),
    panel.background = element_rect(fill = "white",color = "transparent"),
    panel.border = element_rect(fill = "white",color = "transparent"),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.direction = "horizontal",
    legend.position = c(0.5,0.5))-> p11

##
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

library(grid)
svg("output/Fig3.svg",width = 10,height = 14)
(m1/m2/m3/m4)+plot_annotation(subtitle = "a. Counts of transition events")&theme(plot.subtitle = element_text(size = 12, face = "bold")) -> part1
(p1|p2|p3)/
  (p4|p5|p6)/
  (p7|p8|p9) -> part2
p10 -> part3
# viewport(x = 0, y = 0.78, width = 1, height = 0.22, just = c("left", "bottom")) -> vp1
# viewport(x = 0, y = 0, width = 1, height = 0.8, just = c("left", "bottom")) -> vp2
# viewport(x = 0.35, y = 0.01, width = 0.6, height = 0.2, just = c("left", "bottom")) -> vp3
viewport(x = 0, y = 0.72, width = 1, height = 0.28, just = c("left", "bottom")) -> vp1
viewport(x = 0, y = 0.149, width = 1, height = 0.575, just = c("left", "bottom")) -> vp2
viewport(x = 0.008, y = 0, width = 0.425, height = 0.234, just = c("left", "bottom")) -> vp3
viewport(x = 0.55, y = 0.015, width = 0.4, height = 0.1, just = c("left", "bottom")) -> vp4
print(part2,vp = vp2)
print(part1,vp = vp1)
print(p10, vp = vp3)
print(p11, vp = vp4)
dev.off()

pdf("output/Fig3.pdf",width = 10,height = 14)
(m1/m2/m3/m4)+plot_annotation(subtitle = "a. Counts of transition events")&theme(plot.subtitle = element_text(size = 12, face = "bold")) -> part1
(p1|p2|p3)/
  (p4|p5|p6)/
  (p7|p8|p9) -> part2
p10 -> part3
# viewport(x = 0, y = 0.78, width = 1, height = 0.22, just = c("left", "bottom")) -> vp1
# viewport(x = 0, y = 0, width = 1, height = 0.8, just = c("left", "bottom")) -> vp2
# viewport(x = 0.35, y = 0.01, width = 0.6, height = 0.2, just = c("left", "bottom")) -> vp3
viewport(x = 0, y = 0.72, width = 1, height = 0.28, just = c("left", "bottom")) -> vp1
viewport(x = 0, y = 0.149, width = 1, height = 0.575, just = c("left", "bottom")) -> vp2
viewport(x = 0.008, y = 0, width = 0.425, height = 0.234, just = c("left", "bottom")) -> vp3
viewport(x = 0.55, y = 0.015, width = 0.4, height = 0.1, just = c("left", "bottom")) -> vp4
print(part2,vp = vp2)
print(part1,vp = vp1)
print(p10, vp = vp3)
print(p11, vp = vp4)
dev.off()

tiff("output/Fig3.tif",width = 10,height = 14, units = "in", res = 400, compression = "lzw")
(m1/m2/m3/m4)+plot_annotation(subtitle = "a. Counts of transition events")&theme(plot.subtitle = element_text(size = 12, face = "bold")) -> part1
(p1|p2|p3)/
  (p4|p5|p6)/
  (p7|p8|p9) -> part2
p10 -> part3
# viewport(x = 0, y = 0.78, width = 1, height = 0.22, just = c("left", "bottom")) -> vp1
# viewport(x = 0, y = 0, width = 1, height = 0.8, just = c("left", "bottom")) -> vp2
# viewport(x = 0.35, y = 0.01, width = 0.6, height = 0.2, just = c("left", "bottom")) -> vp3
viewport(x = 0, y = 0.72, width = 1, height = 0.28, just = c("left", "bottom")) -> vp1
viewport(x = 0, y = 0.149, width = 1, height = 0.575, just = c("left", "bottom")) -> vp2
viewport(x = 0.008, y = 0, width = 0.425, height = 0.234, just = c("left", "bottom")) -> vp3
viewport(x = 0.55, y = 0.015, width = 0.4, height = 0.1, just = c("left", "bottom")) -> vp4
print(part2,vp = vp2)
print(part1,vp = vp1)
print(p10, vp = vp3)
print(p11, vp = vp4)
dev.off()
