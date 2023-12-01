#==load packages==
library(stringr)
library(Biostrings)
library(ggplot2)
library(tidyverse)
library(readxl)
library(lubridate)
library(ggpubr)
library(rworldmap)
library(ape)
library(ggtree)
library(phangorn)
library(phytools)
library(treeio)
library(zoo)
library(patchwork)
library(ggsci)
library(scales)
library(sf)
library(rgdal)
library(patchwork)
library(grid)

#==define date==
date_match <- data.frame(date = seq(as.Date("1995-01-02"),as.Date("2023-12-25"),7))
date_match$ISO_YEAR <- isoyear(date_match$date)
date_match$ISO_WEEK <- isoweek(date_match$date)

#================================================================
#==Genetic data were not provided due to the policy restriction==
#================================================================
#==influenza A viruses==
#1.read metadata
#GISAID after the de-duplication
gisaid_h1n1_meta <- read.csv("../../../flu_data/h1n1/h1n1_ha_meta_high.csv")
gisaid_h1n1_meta <- gisaid_h1n1_meta[,-12]
names(gisaid_h1n1_meta)[c(7,24)] <- c("alpha.3","seqName")
gisaid_h3n2_meta <- read.csv("../../../flu_data/h3n2/h3n2_ha_meta_high.csv")
gisaid_h3n2_meta <- gisaid_h3n2_meta[,-12]
names(gisaid_h3n2_meta)[c(7,24)] <- c("alpha.3","seqName")
#NCBI after the de-duplication
ncbi_h1n1_meta <- read.csv("../../../flu_data/h1n1/ncbi_ha_meta_h1n1_high.csv")
ncbi_h3n2_meta <- read.csv("../../../flu_data/h3n2/ncbi_ha_meta_h3n2_high.csv")

#2.read sequence
gisaid_h1n1_seq <- readDNAStringSet("../../../flu_data/h1n1/h1n1.aligned_high_qua.fasta")
gisaid_h3n2_seq <- readDNAStringSet("../../../flu_data/h3n2/h3n2.aligned_high_qua.fasta")
ncbi_h1n1_seq <- readDNAStringSet("../../../flu_data/h1n1/ncbi_seq_h1n1_ha_align_high.fasta")
ncbi_h3n2_seq <- readDNAStringSet("../../../flu_data/h3n2/ncbi_seq_h3n2_ha_align_high.fasta")

#3.combine metadata and sequence
h1n1_meta <- rbind(gisaid_h1n1_meta, ncbi_h1n1_meta[,which(names(ncbi_h1n1_meta) %in% names(gisaid_h1n1_meta))]); remove(gisaid_h1n1_meta,ncbi_h1n1_meta)
h3n2_meta <- rbind(gisaid_h3n2_meta, ncbi_h3n2_meta[,which(names(ncbi_h3n2_meta) %in% names(gisaid_h3n2_meta))]); remove(gisaid_h3n2_meta,ncbi_h3n2_meta)
h1n1_meta$source[str_detect(h1n1_meta$Isolate_Id, "EPI")] <- "GISAID"
h1n1_meta$source[!str_detect(h1n1_meta$Isolate_Id, "EPI")] <- "GenBank"
h3n2_meta$source[str_detect(h3n2_meta$Isolate_Id, "EPI")] <- "GISAID"
h3n2_meta$source[!str_detect(h3n2_meta$Isolate_Id, "EPI")] <- "GenBank"
# table(h3n2_meta$source)
# table(h1n1_meta$source)
h1n1_seq <- c(gisaid_h1n1_seq, ncbi_h1n1_seq); remove(gisaid_h1n1_seq,ncbi_h1n1_seq)
h3n2_seq <- c(gisaid_h3n2_seq, ncbi_h3n2_seq); remove(gisaid_h3n2_seq,ncbi_h3n2_seq)
# table(duplicated(h1n1_seq))
# table(duplicated(h3n2_seq))

#==influenza B viruses==
#1.read metadata
#GISAID after the de-duplication
gisaid_bv_meta <- read.csv("../../../flu_data/bv_ha_meta_high.csv")
gisaid_bv_meta <- gisaid_bv_meta[,-12]
names(gisaid_bv_meta)[c(7,23)] <- c("alpha.3","seqName")
gisaid_by_meta <- read.csv("../../../flu_data/by_ha_meta_high.csv")
gisaid_by_meta <- gisaid_by_meta[,-12]
names(gisaid_by_meta)[c(7,23)] <- c("alpha.3","seqName")
#NCBI after the de-duplication
ncbi_bv_meta <- read.csv("../../../flu_data/ncbi_ha_meta_bv_high.csv") %>% select(names(gisaid_bv_meta))
ncbi_by_meta <- read.csv("../../../flu_data/ncbi_ha_meta_by_high.csv") %>% select(names(gisaid_by_meta))

#2.read sequence
gisaid_bv_seq <- readDNAStringSet("../../../flu_data/bv_ha/nextclade.aligned_high_qua.fasta")
gisaid_by_seq <- readDNAStringSet("../../../flu_data/by_ha/nextclade.aligned_high_qua.fasta")
ncbi_bv_seq <- readDNAStringSet("../../../flu_data/ncbi_seq_bv_ha_align_high.fasta")
ncbi_by_seq <- readDNAStringSet("../../../flu_data/ncbi_seq_by_ha_align_high.fasta")

#3.combine metadata and sequence
bv_meta <- rbind(gisaid_bv_meta, ncbi_bv_meta); remove(gisaid_bv_meta,ncbi_bv_meta)
by_meta <- rbind(gisaid_by_meta, ncbi_by_meta); remove(gisaid_by_meta,ncbi_by_meta)
bv_meta$source[str_detect(bv_meta$Isolate_Id, "EPI")] <- "GISAID"
bv_meta$source[!str_detect(bv_meta$Isolate_Id, "EPI")] <- "GenBank"
by_meta$source[str_detect(by_meta$Isolate_Id, "EPI")] <- "GISAID"
by_meta$source[!str_detect(by_meta$Isolate_Id, "EPI")] <- "GenBank"
bv_seq <- c(gisaid_bv_seq, ncbi_bv_seq); remove(gisaid_bv_seq,ncbi_bv_seq)
by_seq <- c(gisaid_by_seq, ncbi_by_seq); remove(gisaid_by_seq,ncbi_by_seq)
# table(duplicated(bv_seq))
# table(duplicated(by_seq))

#==Epidemiological surveillance data==
epi_final <- read.csv("../data_part/epi_data/epi_final_no_roll.csv") %>% mutate(date = as.Date(date))
epi_glo_flu_region1 <- read.csv("../data_part/epi_data/epi_glo_flu_region1_no_roll.csv") %>% mutate(date = as.Date(date))
epi_glo_flua1 <- read.csv("../data_part/epi_data/epi_glo_flua1_no_roll.csv") %>% mutate(date = as.Date(date))
epi_glo_flub1 <- read.csv("../data_part/epi_data/epi_glo_flub1_no_roll.csv") %>% mutate(date = as.Date(date))

#========
#==plot==
#========
h1n1_meta$date <- as.Date(h1n1_meta$date)
h3n2_meta$date <- as.Date(h3n2_meta$date)
bv_meta$date <- as.Date(bv_meta$date)
by_meta$date <- as.Date(by_meta$date)

h1n1_num <- h1n1_meta %>%
  group_by(date) %>%
  summarise(h1n1_seq = n())
h3n2_num <- h3n2_meta %>%
  group_by(date) %>%
  summarise(h3n2_seq = n())
bv_num <- bv_meta %>%
  group_by(date) %>%
  summarise(bv_seq = n())
by_num <- by_meta %>%
  group_by(date) %>%
  summarise(by_seq = n())

epi_glo_flua1 <- left_join(epi_glo_flua1, h1n1_num) %>% left_join(h3n2_num)
epi_glo_flub1 <- left_join(epi_glo_flub1, bv_num) %>% left_join(by_num)

#==Figure 1==
colors <- c(pal_npg("nrc", alpha =1)(10)[c(1:7,9:10)],"darkred","#FADDA9","grey80")
colors1 <- c(pal_aaas("default", alpha =0.7)(10))
show_col(colors)
show_col(colors1)
value = c("Japan/Korea" = colors[1],"Western Asia" = colors[3],"Northern America" = colors[6],
            "South-eastern Asia"= colors[4], "Southern Asia"= colors[5], "Europe"= colors[2], "Oceania"= colors[7],
            "North China" = colors[8], "South China" = colors[11], "Russia"= colors[10],  "Southern America"= colors[12], 
           "Africa"= colors[9], "Americas" = colors1[7], "Asia" = "#D8BFD8", "China" = colors1[4])

max(epi_glo_flua1$h1n1_num[epi_glo_flua1$date >= as.Date("2020-05-01") & epi_glo_flua1$date <= as.Date("2021-07-31")]/
  epi_glo_flua1$test[epi_glo_flua1$date >= as.Date("2020-05-01") & epi_glo_flua1$date <= as.Date("2021-07-31")])
max(epi_glo_flua1$h3n2_num[epi_glo_flua1$date >= as.Date("2020-05-01") & epi_glo_flua1$date <= as.Date("2021-07-31")]/
      epi_glo_flua1$test[epi_glo_flua1$date >= as.Date("2020-05-01") & epi_glo_flua1$date <= as.Date("2021-07-31")])
max(epi_glo_flub1$BV_num[epi_glo_flua1$date >= as.Date("2020-05-01") & epi_glo_flua1$date <= as.Date("2021-07-31")]/
      epi_glo_flub1$test[epi_glo_flua1$date >= as.Date("2020-05-01") & epi_glo_flua1$date <= as.Date("2021-07-31")])
max(epi_glo_flub1$BY_num[epi_glo_flua1$date >= as.Date("2020-05-01") & epi_glo_flua1$date <= as.Date("2021-07-31")]/
      epi_glo_flub1$test[epi_glo_flua1$date >= as.Date("2020-05-01") & epi_glo_flua1$date <= as.Date("2021-07-31")])

ggplot(data = h1n1_meta) + 
  annotate("rect", xmin = as.Date("2020-02-01"),xmax = as.Date("2021-07-31"),
           ymin = 0,ymax = 1300,alpha = 0.2,fill = colors[5])+
  geom_bar(aes(x = date, fill = region))+
  geom_ribbon(data = epi_glo_flua1, aes(x = date,
                                        ymin = h1n1_num_roll_LL/test_roll*6000,
                                        ymax = h1n1_num_roll_UL/test_roll*6000),
              stat = "identity",fill = "darkred", alpha = 0.4)+
  geom_line(data = epi_glo_flua1, aes(x = date, y = h1n1_num_roll/test_roll*6000),
            stat = "identity",size = 0.5,color = "darkred")+
  # annotate("text", x= c(as.Date("2020-11-01")),y = 750, size = 3,
  #          label = c("Pandemic\nperiod"))+
  scale_x_date("", date_breaks = "1 year",date_labels = "%Y",
               limits = c(as.Date("2015-10-01"),as.Date("2023-08-31")))+
  scale_y_continuous("No. of genomes", 
                     limits = c(-20,1300),
                     breaks = seq(0,1300,300),
                     expand = c(0,0),
                     sec.axis = sec_axis(trans=~.* 1/6000,
                                         name="Positivity rate (%)",
                                         breaks = seq(0,0.2,0.05),
                                         labels = seq(0,20,5)
                     ))+
  theme_bw()+
  scale_fill_manual("Continents",values = value)+
  theme(legend.position = c(0.75,0.88))+
  theme(axis.text.y.right = element_text(color = "darkred"),
        panel.grid.major.x  = element_blank(),
        panel.grid.minor.x  = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.3,"cm"),
        legend.background = element_blank(),
        axis.title.y.right = element_text(color = "darkred"),
        axis.title.x = element_blank())+
  guides(fill = guide_legend(nrow = 1))+
  labs(subtitle = "H1N1pdm09") -> p1_2

ggplot(data = h3n2_meta) + 
  annotate("rect", xmin = as.Date("2020-02-01"),xmax = as.Date("2021-07-31"),
           ymin = 0,ymax = 1900,alpha = 0.2,fill = colors[5])+
  geom_bar(aes(x = date, fill = region))+
  geom_ribbon(data = epi_glo_flua1, aes(x = date, 
                                        ymin = h3n2_num_roll_LL/test_roll*8000,
                                        ymax = h3n2_num_roll_UL/test_roll*8000),
              stat = "identity",fill = "darkred", alpha = 0.4)+
  geom_line(data = epi_glo_flua1 ,
            aes(x = date, y = h3n2_num_roll/test_roll*8000), 
            stat = "identity",size = 0.5,color = "darkred")+
  scale_x_date("", date_breaks = "1 year",date_labels = "%Y",
               limits = c(as.Date("2015-10-01"),as.Date("2023-08-30")))+
  scale_y_continuous("No. of genomes", 
                     limits = c(-40,1900),
                     breaks = seq(0,1600,400),
                     expand = c(0,0),
                     sec.axis = sec_axis(trans=~.* 1/8000,
                                         name="Positivity rate (%)",
                                         breaks = seq(0,0.2,0.05),
                                         labels = seq(0,20,5)
                     ))+
  theme_bw()+
  # scale_fill_manual("Genome sources",values = c("darkorange","lightblue"))+
  scale_fill_manual("Continents",values = value)+
  theme(legend.position = 'none')+
  theme(axis.text.y.right = element_text(color = "darkred"),
        axis.title.x = element_blank(),
        panel.grid.major.x  = element_blank(),
        panel.grid.minor.x  = element_blank(),
        axis.title.y.right = element_text(color = "darkred"))+
  labs(subtitle = "H3N2")-> p2_2

bv_meta$date <-as.Date(bv_meta$date)
ggplot(data = bv_meta) + 
  # annotate("rect", xmin = as.Date("2019-05-01"),xmax = as.Date("2020-02-28"),
  #          ymin = 0,ymax = 450,alpha = 0.1)+
  annotate("rect", xmin = as.Date("2020-02-01"),xmax = as.Date("2021-07-31"),
           ymin = 0,ymax = 600,alpha = 0.2,fill = colors[5])+
  # annotate("rect", xmin = as.Date("2022-09-01"),xmax = as.Date("2023-06-30"),
  #          ymin = 0,ymax = 450,alpha = 0.1,fill = "green")+
  geom_bar(aes(x = date, fill = region))+
  geom_ribbon(data = epi_glo_flub1, aes(x = date, ymin = BV_num_roll_LL/test_roll*5000,
                                        ymax = BV_num_roll_UL/test_roll*5000),
              stat = "identity",fill = "darkred", alpha = 0.4)+
  geom_line(data = epi_glo_flub1 ,
            aes(x = date, y = BV_num_roll/test_roll*5000), stat = "identity",size = 0.5,color = "darkred")+
  # geom_line(data = epi_glo_flub1 ,
  #           aes(x = date, y = BV_num_roll/30), stat = "identity",size = 0.2,color = "blue")+
  scale_x_date("", date_breaks = "1 year",date_labels = "%Y",
               limits = c(as.Date("2015-10-01"),as.Date("2023-08-30")))+
  scale_y_continuous("No. of genomes", 
                     breaks = seq(0,600,150),
                     limits = c(-10,600),
                    
                     expand = c(0,0),
                     sec.axis = sec_axis(trans=~.* 1/5000,
                                         name="Positivity rate (%)",
                                         breaks = seq(0,0.12,0.03),
                                         labels = seq(0,12,3)
                     ))+
  theme_bw()+
  # scale_fill_manual("Genome sources",values = c("darkorange","lightblue"))+
  scale_fill_manual("Continents",values = value)+
  theme(axis.text.y.right = element_text(color = "darkred"),
        axis.title.x = element_blank(),
        panel.grid.major.x  = element_blank(),
        panel.grid.minor.x  = element_blank(),
        axis.title.y.right = element_text(color = "darkred"))+
  theme(legend.position = "none")+
  labs(subtitle = "B/Victoria")-> p3_2

by_meta$date <-as.Date(by_meta$date)
ggplot(data = by_meta) + 
  # annotate("rect", xmin = as.Date("2019-05-01"),xmax = as.Date("2020-02-28"),
  #          ymin = 0,ymax = 600,alpha = 0.1)+
  annotate("rect", xmin = as.Date("2020-02-01"),xmax = as.Date("2021-07-31"),
           ymin = 0,ymax = 610,alpha = 0.2,fill = colors[5])+
  # annotate("rect", xmin = as.Date("2022-09-01"),xmax = as.Date("2023-06-30"),
  #          ymin = 0,ymax = 600,alpha = 0.1,fill = "green")+
  geom_bar(aes(x = date, fill = region))+
  geom_ribbon(data = epi_glo_flub1[epi_glo_flub1$date < as.Date("2020-07-01"),],
              aes(x = date, ymin = BY_num_roll_LL/test_roll*4000,
                                        ymax = BY_num_roll_UL/test_roll*4000),
              stat = "identity",fill = "darkred", alpha = 0.4)+
  geom_line(data = epi_glo_flub1[epi_glo_flub1$date < as.Date("2020-07-01"),] ,
            aes(x = date, y = BY_num_roll/test_roll*4000), stat = "identity",size = 0.5,color = "darkred")+
  # geom_line(data = epi_glo_flub1 ,
  #           aes(x = date, y = BY_num_roll/40), stat = "identity",size = 0.2,color = "blue")+
  scale_x_date("Date of collection", date_breaks = "1 year",date_labels = "%Y",
               limits = c(as.Date("2015-10-01"),as.Date("2023-08-30")))+
  scale_y_continuous("No. of genomes", 
                     breaks = seq(0,600,200),
                     limits = c(-15,610),
                      expand = c(0,0),
                     sec.axis = sec_axis(trans=~.* 1/4000,
                                         name="Positivity rate (%)",
                                         breaks = seq(0,0.15,0.05),
                                         labels = seq(0,15,5)
                     ))+
  theme_bw()+
  # scale_fill_manual("Genome sources",values = c("darkorange","lightblue"))+
  scale_fill_manual("Continents",values = value)+
  theme(legend.position = "none")+
  theme(axis.text.y.right = element_text(color = "darkred"),
        panel.grid.major.x  = element_blank(),
        panel.grid.minor.x  = element_blank(),
        axis.title.y.right = element_text(color = "darkred"))+
  labs(subtitle = "B/Yamagata")-> p4_2

#==two surveillance indicators==
ggplot(data = epi_glo_flua1[epi_glo_flua1$date >= as.Date("2016-08-01"),]) +
  annotate("rect", xmin = as.Date("2020-02-01"),xmax = as.Date("2021-07-31"),
           ymin = 0,ymax = 400000,alpha = 0.2,fill = colors[5])+
  annotate("text", x= c(as.Date("2020-11-01")),y = 320000, size = 3,
           label = c("High-NPIs phase\nof the pandemic"))+
  geom_line(aes(x = date, y = test_roll))+
  scale_x_date("Date", date_breaks = "1 year", date_labels = "%Y")+
  scale_y_continuous("No. of specimens (× 1k)", labels = seq(0,400,100),
                     limits = c(0, 400000),expand = c(0,0))+
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),)+
  labs(subtitle = "Intensity of virological surveillance")-> fig1

meta <- rbind(h1n1_meta[,c(23,24)], h3n2_meta[,c(23,24)], bv_meta[,c(22, 23)], by_meta[,c(22, 23)])
seq_epi <- meta %>%
  group_by(date) %>%
  summarise(seq_no = n()) %>%
  left_join(epi_glo_flua1) %>%
  left_join(epi_glo_flub1[,c(2:15)]) %>%
  ungroup() %>%
  group_by(date) %>%
  summarise(epi_num = sum(h1n1_num+h3n2_num+BV_num+BY_num, na.rm = T),
            seq_no = sum(seq_no, na.rm = T))

seq_epi1 <- seq_epi %>%
  mutate(epi_num_roll = rollmean(epi_num, k=5, fill=NA, align='center'),
         seq_no_roll = rollmean(seq_no, k=5, fill=NA, align='center'))
ggplot(data = seq_epi1[seq_epi1$date >= as.Date("2016-08-01"),]) +
  annotate("rect", xmin = as.Date("2020-02-01"),xmax = as.Date("2021-07-31"),
           ymin = 0,ymax = 0.3,alpha = 0.2,fill = colors[5])+
  # geom_hline(yintercept = 0.05, linetype = 2, color = "red")+
  geom_line(aes(x = date, y = seq_no_roll/epi_num_roll))+
  scale_x_date("Date", date_breaks = "1 year", date_labels = "%Y")+
  scale_y_continuous("Proportion of cases\nsequenced (%)",
                     expand = c(0,0),limits = c(0,0.3),
                     breaks = seq(0,0.3,0.1),labels = seq(0, 30, 10))+
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        panel.background = element_rect(fill = "transparent", color = "transparent"),
        plot.background = element_rect(fill = "transparent", color = "transparent"))+
  labs(subtitle = "Intensity of genomic surveillance") -> fig2

#==plot-by-regions==
epi_glo_flu_region2 <- epi_glo_flu_region1[epi_glo_flu_region1$region_final %in% c("Africa", "China","South-eastern Asia", "Southern Asia"),]
epi_glo_flu_region2 <- epi_glo_flu_region2[(epi_glo_flu_region2$date >= as.Date("2019-10-01") & 
                                             epi_glo_flu_region2$date < as.Date("2021-10-01")),]
ggplot(epi_glo_flu_region2) +
  annotate("rect", xmin = as.Date("2020-02-01"),xmax = as.Date("2021-07-31"),
           ymin = 0,ymax = 0.3,alpha = 0.2,fill = colors[5])+
  # geom_ribbon(aes(x = date, ymin = h1n1_num_roll_LL/test_roll,ymax = h1n1_num_roll_UL/test_roll, color = region_final),
  #             stat = "identity",fill = "darkred", alpha = 0.4) +
  geom_line(aes(x = date, y = h1n1_num_roll/test_roll, color = region_final))+
  scale_x_date("Date", date_breaks = "8 month",date_labels = "%b %Y")+
  scale_y_continuous("Positivity rate (%)",breaks = seq(0,0.3,0.1),
                     labels = seq(0,30,10), expand = c(0, 0),limits = c(-0.01,0.3))+
  theme_bw()+
  theme(
    # axis.title.y = element_blank(),
        # axis.text.y = element_text(color = "darkred")
        )+
  scale_color_manual("Regions",values = value,
                     labels = c("Africa", "China", "Southeast Asia", "South Asia"))+
  theme(axis.title.x = element_blank(),
        panel.grid.major.x  = element_blank(),
        panel.grid.minor.x  = element_blank(),
        # legend.background = element_blank(),
        legend.key = element_rect(fill = "transparent", color ="transparent"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 7),
        legend.key.size = unit(0.3,"cm"),
        # legend.title = element_blank(),
        # axis.title.y.right = element_text(color = "darkred")
        )+
  theme(legend.position = c(0.6,0.65))-> p9

ggplot(epi_glo_flu_region2) +
  annotate("rect", xmin = as.Date("2020-02-01"),xmax = as.Date("2021-07-31"),
           ymin = 0,ymax = 0.3,alpha = 0.2,fill = colors[5])+
  geom_line(aes(x = date, y = h3n2_num_roll/test_roll, color = region_final))+
  scale_x_date("Date", date_breaks = "8 month",date_labels = "%b %Y")+
  scale_y_continuous("Positivity rate (%)",breaks = seq(0,0.3,0.1),
                     labels = seq(0,30,10), expand = c(0, 0),limits = c(-0.01,0.3))+
  theme_bw()+
  theme(
    # axis.title.y = element_blank(),
        # axis.text.y = element_text(color = "darkred")
        )+
  theme(axis.title.x = element_blank(),
        panel.grid.major.x  = element_blank(),
        panel.grid.minor.x  = element_blank(),
        legend.title = element_blank(),
        # axis.title.y.right = element_text(color = "darkred")
        )+
  scale_color_manual("Regions",values = value)+
  theme(legend.position = "none") -> p10

ggplot(epi_glo_flu_region2) +
  annotate("rect", xmin = as.Date("2020-02-01"),xmax = as.Date("2021-07-31"),
           ymin = 0,ymax = 0.15,alpha = 0.2,fill = colors[5])+
  geom_line(aes(x = date, y = BV_num_roll/test_roll, color = region_final))+
  scale_x_date("Date", date_breaks = "8 month",date_labels = "%b %Y")+
  scale_y_continuous("Positivity rate (%)",breaks = seq(0,0.15,0.05),
                     labels = seq(0,15,5), expand = c(0, 0),limits = c(-0.01/2,0.3/2))+
  theme_bw()+
  theme(
    # axis.title.y = element_blank(),
        # axis.text.y = element_text(color = "darkred")
        )+
  theme(axis.title.x = element_blank(),
        panel.grid.major.x  = element_blank(),
        panel.grid.minor.x  = element_blank(),
        legend.title = element_blank(),
        axis.title.y.right = element_text(color = "darkred"))+
  scale_color_manual("Regions",values = value)+
  theme(legend.position = "none")-> p11

ggplot()+theme(panel.background = element_rect(fill = "transparent", color = "transparent"))+
  labs(y = "")-> f_n

#==map plot==
#=========================================================================
#we didn't provide the map data as below due to the data sharing agreement
#=========================================================================
ggplot() +
  geom_sf(data = worldmap0, color = "grey50", fill = "grey92")+
  geom_sf(data = worldmap1[worldmap1$region_final %in% c("Southern Asia","South-eastern Asia","Africa","China"),],aes(fill = region_final), color = "black", size = 0.1)+
  geom_sf(data = nine, color="black",size = 0.01)+
  annotate("rect", xmin = -180, xmax = 180, ymin = -23.27, ymax = 23.27,size = 0.1, alpha = 0.3, fill = "lightblue")+
  geom_hline(yintercept = c(-35, 35), size = 0.1, linetype = 2, color = "red")+
  annotate("text", x = -180, y = 0, label = "Tropics", hjust = 0, size = 2)+
  annotate("text", x = -180, y = -29, label = "Subtropics", hjust = 0, size = 2)+
  annotate("text", x = -180, y = 29, label = "Subtropics", hjust = 0, size = 2)+
  theme_void()+
  theme(legend.position = "top")+
  guides(fill = F)+
  scale_fill_manual("Geographic regions",values = value) -> map

#==output==
pdf("../output/Fig1.pdf",height = 10, width = 10)
(fig1|fig2)/
  ((p1_2|p9)+plot_layout(nrow = 1, ncol = 2,widths = c(1,0.4)))/
  ((p2_2|p10)+plot_layout(nrow = 1, ncol = 2,widths = c(1,0.402)))/
     ((p3_2|p11)+plot_layout(nrow = 1, ncol = 2,widths = c(1,0.4)))/
        ((p4_2|f_n)+plot_layout(nrow = 1, ncol = 2,widths = c(1,0.442)))/
           plot_annotation(tag_levels = "a")
viewport(x = 0.71, y = 0.015, width = 0.3, height = 0.18, just = c("left", "bottom")) -> vp1
print(map,vp = vp1)
dev.off()

tiff("../output/Fig1.tif", width = 10, height = 10,
      units = "in", compression = "lzw", res = 400)
(fig1|fig2)/
  ((p1_2|p9)+plot_layout(nrow = 1, ncol = 2,widths = c(1,0.4)))/
  ((p2_2|p10)+plot_layout(nrow = 1, ncol = 2,widths = c(1,0.402)))/
  ((p3_2|p11)+plot_layout(nrow = 1, ncol = 2,widths = c(1,0.4)))/
  ((p4_2|f_n)+plot_layout(nrow = 1, ncol = 2,widths = c(1,0.442)))/
  plot_annotation(tag_levels = "a")
viewport(x = 0.71, y = 0.015, width = 0.3, height = 0.18, just = c("left", "bottom")) -> vp1
print(map,vp = vp1)
dev.off()

svg("../output/Fig1.svg",height = 10, width = 10)
(fig1|fig2)/
  ((p1_2|p9)+plot_layout(nrow = 1, ncol = 2,widths = c(1,0.4)))/
  ((p2_2|p10)+plot_layout(nrow = 1, ncol = 2,widths = c(1,0.402)))/
  ((p3_2|p11)+plot_layout(nrow = 1, ncol = 2,widths = c(1,0.4)))/
  ((p4_2|f_n)+plot_layout(nrow = 1, ncol = 2,widths = c(1,0.442)))/
  plot_annotation(tag_levels = "a")
viewport(x = 0.71, y = 0.015, width = 0.3, height = 0.18, just = c("left", "bottom")) -> vp1
print(map,vp = vp1)
dev.off()