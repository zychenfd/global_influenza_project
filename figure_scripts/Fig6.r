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
library(castor)
library(lubridate)
library(patchwork)
library(grid)

#==pathway and color==
setwd("C:/Users/zyche/Nutstore/1/Evolution_study/Flu_phylogeography/analyses/scripts")
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
          "Africa"= colors[9], "Americas" = colors1[1], "Asia" = colors1[2], "China" = colors1[4])

#==1. effective pop size
h1n1_even <- read.delim("../phylogeography/4.1post-analyses/pop_size/h1n1_even_pop_size.txt") %>%
  mutate(date = as.Date(date), type = "H1N1pdm09") %>%
  filter(date >= as.Date("2012-01-01"))
h3n2_even <- read.delim("../phylogeography/4.1post-analyses/pop_size/h3n2_even_pop_size.txt") %>%
  mutate(date = as.Date(date), type = "H3N2") %>%
  filter(date >= as.Date("2012-01-01"))
bv_even <- read.delim("../phylogeography/4.1post-analyses/pop_size/bv_even_pop_size.txt") %>%
  mutate(date = as.Date(date), type = "B/Victoria") %>%
  filter(date >= as.Date("2012-01-01"))
by_even <- read.delim("../phylogeography/4.1post-analyses/pop_size/by_even_pop_size.txt") %>%
  mutate(date = as.Date(date), type = "B/Yamagata") %>%
  filter(date >= as.Date("2012-01-01"))
pop_size <- rbind(h1n1_even, h3n2_even, bv_even, by_even)
factor(pop_size$type, levels = unique(pop_size$type)) -> pop_size$type

ggplot(data = pop_size) +
  annotate("rect", xmin = as.Date("2020-02-01"),xmax = as.Date("2021-07-31"),
           ymin = 0,ymax = 300,alpha = 0.2,fill = colors[5])+
  # geom_vline(xintercept = c(as.Date("2020-02-01"),as.Date("2021-08-01")), linetype = 2)+
  geom_ribbon(aes(x = date, ymin = lower, ymax = upper, fill = type), alpha = 0.2)+
  geom_line(aes(x = date, y = mean, color = type, group = type))+
  geom_point(aes(x = date, y = mean, fill = type), shape = 21)+
  scale_y_log10(expand = c(0,0))+
  # scale_y_continuous(expand = c(0,0))+
  scale_x_date(expand = c(0,0),date_labels = "%Y",date_breaks = "1 year",
               limits = c(as.Date("2012-07-01"),as.Date("2023-07-30")))+
  # coord_cartesian(ylim = c(0,151))+
  theme_bw()+
  theme(legend.position = "top",
        axis.title.x = element_blank(),
        legend.background = element_blank())+
  labs(x = "Year", y = "Relative genetic diversity", tag = "b")+
  # scale_color_manual("",values = colors1[c(1,5)],
  #                    labels = c("B/Victoria (V1A.3 sub-clades)","B/Yamagata (Y3 clades with M251V mutation)"))+
  # scale_fill_manual("",values = colors1[c(1,5)],
  #                   labels = c("B/Victoria (V1A.3 sub-clades)","B/Yamagata (Y3 clades with M251V mutation)"))+
  scale_color_manual("",values = colors1[c(1:3,9)])+
  scale_fill_manual("",values = colors1[c(1:3,9)])-> p2

#==2. MCC tree
#B/Yamagata
tree_by <- read.beast("../phylogeography/4.1post-analyses/mcc_tree/by_from2011_mcc.tre")
tree_by@phylo$edge.length[tree_by@phylo$edge.length < 0] <- 0.00001
by_even_meta <- read.csv("../../../flu_data/Final_data1/by_tree_files/metadata_by_even1.csv") %>% select(c("seqName", "region_final"))
names(by_even_meta)[1] <- "taxa"
ggtree(tree_by, mrsd = as.Date("2020-03-20"), as.Date=TRUE,color='grey40',size=0.1) %<+% by_even_meta + geom_tippoint(aes(fill = region_final),size=2.5, color='black',shape=21, stroke=0.1)+
  # geom_vline(xintercept = as.Date("2020-02-01"), linetype = 2, size = 0.5, color = "red")+
  annotate("rect", xmin = as.Date("2020-02-01"),xmax = as.Date("2021-07-31"),
           ymin = -50,ymax = 2050,alpha = 0.2,fill = colors[5])+
  theme_tree2() +
  scale_x_date("Date",date_labels = "%Y",date_breaks = "2 year", expand = c(0.02,0), limits = c(as.Date("2006-01-01"), as.Date("2023-06-30")))+
  scale_fill_manual("", values = value)+
  theme(legend.position = "none")+
  scale_y_continuous(expand = c(0,0), limits = c(-50, 2050))+
  theme(
    # axis.text.x =  element_text(angle = 45,hjust =1,vjust = 1),
    axis.title.x = element_blank(),
    panel.border = element_rect(fill = "transparent", color = "transparent"),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid.major.x = element_line(color = "grey90"))+
  # annotate("text", x = as.Date("2007-12-01"), y = 2000, label = "B/Yamagata",  size = 4)
  labs(subtitle  = "B/Yamagata", tag = "a")-> p1

by_even_meta1 <- read.csv("../../../flu_data/Final_data1/by_tree_files/metadata_by_even1.csv") %>% select(c("seqName", "clade"))
names(by_even_meta1)[1] <- "taxa"
ggtree(tree_by, mrsd = as.Date("2020-03-20"), as.Date=TRUE,color='grey40',size=0.1) %<+% 
  by_even_meta1 + geom_tippoint(aes(fill = clade),size=1.5, color='black',shape=21, stroke=0.1)+
  # geom_vline(xintercept = as.Date("2020-02-01"), linetype = 2, size = 0.5, color = "red")+
  theme_tree2() +
  scale_x_date("Date",date_labels = "%Y",date_breaks = "3 year", expand = c(0.02,0))+
  scale_fill_manual("Clades", values = colors1[c(2,5,4)])+
  theme(legend.position = c(0.2,0.7))+
  scale_y_continuous(expand = c(0.02,0))+
  theme(
    # axis.text.x =  element_text(angle = 45,hjust =1,vjust = 1),
    axis.title.x = element_blank(),
    panel.border = element_rect(fill = "transparent", color = "transparent"),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()) -> p1_1

ggtree(tree_by, mrsd = as.Date("2020-03-20"), as.Date=TRUE,color='grey40',size=0.1) %<+% by_even_meta + geom_tippoint(aes(fill = region_final),size=2.5, color='black',shape=21, stroke=0.1)+
  # geom_vline(xintercept = as.Date("2020-02-01"), linetype = 2, size = 0.5, color = "red")+
  annotate("rect", xmin = as.Date("2020-02-01"),xmax = as.Date("2021-07-31"),
           ymin = -50,ymax = 2050,alpha = 0.2,fill = colors[5])+
  theme_tree2() +
  scale_x_date("Date",date_labels = "%Y",date_breaks = "2 year", expand = c(0.02,0), limits = c(as.Date("2006-01-01"), as.Date("2023-06-30")))+
  scale_fill_manual("", values = value)+
  theme(legend.position = "none")+
  scale_y_continuous(expand = c(0,0), limits = c(-50, 2050))+
  theme(
    # axis.text.x =  element_text(angle = 45,hjust =1,vjust = 1),
    axis.title.x = element_blank(),
    panel.border = element_rect(fill = "transparent", color = "transparent"),
    # plot.subtitle = element_text(hjust = 0.5),
    panel.grid.major.x = element_line(color = "grey90"))+
  # annotate("text", x = as.Date("2007-12-01"), y = 2000, label = "B/Yamagata",  size = 4)
  labs(subtitle  = "B/Yamagata")-> p1_2

#H3N2
tree_h3n2 <- read.beast("../phylogeography/4.1post-analyses/mcc_tree/h3n2_from2011_mcc.tre")
tree_h3n2@phylo$edge.length[tree_h3n2@phylo$edge.length < 0] <- 0.00001
h3n2_even_meta <- read.csv("../../../flu_data/Final_data1/h3n2_tree_files/metadata_h3n2_even1.csv") %>% select(c("seqName", "region_final"))
names(h3n2_even_meta)[1] <- "taxa"
ggtree(tree_h3n2, mrsd = as.Date("2023-07-13"), as.Date=TRUE,color='grey40',size=0.1) %<+% h3n2_even_meta + geom_tippoint(aes(fill = region_final),size=2.5, color='black',shape=21, stroke=0.1)+
  # geom_vline(xintercept = as.Date("2020-02-01"), linetype = 2, size = 0.5, color = "red")+
  annotate("rect", xmin = as.Date("2020-02-01"),xmax = as.Date("2021-07-31"),
           ymin = -50,ymax = 2250,alpha = 0.2,fill = colors[5])+
  theme_tree2() +
  scale_x_date("Date",date_labels = "%Y",date_breaks = "2 year", expand = c(0.02,0))+
  scale_fill_manual("", values = value)+
  theme(legend.position = "none")+
  scale_y_continuous(expand = c(0,0), limits = c(-50, 2250))+
  theme(
    # axis.text.x =  element_text(angle = 45,hjust =1,vjust = 1),
    axis.title.x = element_blank(),
    panel.border = element_rect(fill = "transparent", color = "transparent"),
    panel.grid.major.x = element_line(color = "grey90"))+
  # annotate("text", x = as.Date("2007-12-01"), y = 2000, label = "B/Yamagata",  size = 4)
  labs(subtitle  = "H3N2")-> p1_h3n2

#H1N1pdm09
tree_h1n1 <- read.beast("../phylogeography/4.1post-analyses/mcc_tree/h1n1_from2011_mcc.tre")
tree_h1n1@phylo$edge.length[tree_h1n1@phylo$edge.length < 0] <- 0.00001
h1n1_even_meta <- read.csv("../../../flu_data/Final_data1/h1n1_tree_files/metadata_h1n1_even1.csv") %>% select(c("seqName", "region_final"))
names(h1n1_even_meta)[1] <- "taxa"
ggtree(tree_h1n1, mrsd = as.Date("2023-07-26"), as.Date=TRUE,color='grey40',size=0.1) %<+% h1n1_even_meta + geom_tippoint(aes(fill = region_final),size=2.5, color='black',shape=21, stroke=0.1)+
  # geom_vline(xintercept = as.Date("2020-02-01"), linetype = 2, size = 0.5, color = "red")+
  annotate("rect", xmin = as.Date("2020-02-01"),xmax = as.Date("2021-07-31"),
           ymin = -50,ymax = 2100,alpha = 0.2,fill = colors[5])+
  theme_tree2() +
  scale_x_date("Date",date_labels = "%Y",date_breaks = "2 year", expand = c(0.02,0))+
  scale_fill_manual("", values = value)+
  theme(legend.position = "none")+
  scale_y_continuous(expand = c(0,0), limits = c(-50, 2100))+
  theme(
    # axis.text.x =  element_text(angle = 45,hjust =1,vjust = 1),
    axis.title.x = element_blank(),
    panel.border = element_rect(fill = "transparent", color = "transparent"),
    panel.grid.major.x = element_line(color = "grey90"))+
  # annotate("text", x = as.Date("2007-12-01"), y = 2000, label = "B/Yamagata",  size = 4)
  labs(subtitle  = "H1N1pdm09")-> p1_h1n1

#B/Victoria
tree_bv <- read.beast("../phylogeography/4.1post-analyses/mcc_tree/bv_from2011_mcc.tre")
tree_bv@phylo$edge.length[tree_bv@phylo$edge.length < 0] <- 0.00001
bv_even_meta <- read.csv("../../../flu_data/Final_data1/bv_tree_files/metadata_bv_even1.csv") %>% select(c("seqName", "region_final"))
names(bv_even_meta)[1] <- "taxa"
ggtree(tree_bv, mrsd = as.Date("2023-07-24"), as.Date=TRUE,color='grey40',size=0.1) %<+% bv_even_meta + geom_tippoint(aes(fill = region_final),size=2.5, color='black',shape=21, stroke=0.1)+
  # geom_vline(xintercept = as.Date("2020-02-01"), linetype = 2, size = 0.5, color = "red")+
  annotate("rect", xmin = as.Date("2020-02-01"),xmax = as.Date("2021-07-31"),
           ymin = -50,ymax = 2200,alpha = 0.2,fill = colors[5])+
  theme_tree2() +
  scale_x_date("Date",date_labels = "%Y",date_breaks = "2 year", expand = c(0.02,0))+
  scale_fill_manual("", values = value)+
  theme(legend.position = "none")+
  scale_y_continuous(expand = c(0,0), limits = c(-50, 2200))+
  theme(
    # axis.text.x =  element_text(angle = 45,hjust =1,vjust = 1),
    axis.title.x = element_blank(),
    panel.border = element_rect(fill = "transparent", color = "transparent"),
    panel.grid.major.x = element_line(color = "grey90"))+
  # annotate("text", x = as.Date("2007-12-01"), y = 2000, label = "B/Yamagata",  size = 4)
  labs(subtitle  = "B/Victoria")-> p1_bv

#==output for MCC tree (Extended Data figure)==
pdf("output/MCC_tree.pdf", width = 10, height = 10)
(p1_h1n1/p1_bv)|(p1_h3n2/p1_2)
dev.off()

svg("output/MCC_tree.svg", width = 10, height = 10)
(p1_h1n1/p1_bv)|(p1_h3n2/p1_2)
dev.off()

#==3.lineage turnover
by_even1 <- read.csv("../../../flu_data/Final_data1/by_tree_files/metadata_by_even1_hemi.csv") %>% mutate(Collection_Date1 = as.Date(Collection_Date1))
max(as.Date(by_even1$Collection_Date1))
cut(by_even1$Collection_Date1, breaks = c(as.Date("2010-06-30"),as.Date("2011-06-30"),as.Date("2012-06-30"),
                                          as.Date("2013-06-30"),as.Date("2014-06-30"),as.Date("2015-06-30"),
                                          as.Date("2016-06-30"),as.Date("2017-06-30"),
                                          as.Date("2018-06-30"),as.Date("2019-06-30"),as.Date("2020-06-30"),
                                          as.Date("2021-06-30"),as.Date("2022-06-30"),as.Date("2023-06-30"),as.Date("2024-06-30")),
    labels = c(2011:2024), right = T) -> by_even1$North_hem_season
by_even1$South_hem_season <- as.numeric(paste0(substr(by_even1$Collection_Date1,1,4),".5"))
by_even1$season_hemi <- paste0(by_even1$type,"_",by_even1$North_hem_season)
by_even1$season_hemi[by_even1$type == "S_Hemisphere"] <- paste0(by_even1$type[by_even1$type == "S_Hemisphere"],"_",by_even1$South_hem_season[by_even1$type == "S_Hemisphere"])
table(by_even1$season_hemi)

by_trees <- read.beast("../phylogeography/1.3empirical_trees/by_from2011_cut.trees") 
mrac_data1 <- c()
mrac_data1 <- as.data.frame(mrac_data1)
for (i in 1:length(by_trees)) {
  by_trees_1 <- by_trees[[i]]
  ggtree(by_trees_1,mrsd = as.Date("2020-03-20"))+theme_tree2() -> p_tmp
  data <- p_tmp$data
  
  tip_label <- unique(by_even1$season_hemi)
  
  for (j in 1:length(tip_label)) {
    mrca_node <- get_mrca_of_set(by_trees_1@phylo, by_even1$seqName[by_even1$season_hemi == tip_label[j]])
    mrca_time <- data$x[data$node == mrca_node]
    mrac_data <- data.frame(virus = "by", tree_id = i, season = tip_label[j], tmrca = mrca_time)
    mrac_data1 <- rbind(mrac_data1, mrac_data)
  }
}
# write.csv(mrac_data1,"output/by_tmrca_hemi1.csv", row.names = F)
# by_tmrca_hemi <- read.csv("output/by_tmrca_hemi1.csv")

summarise_hpd_lower <- function(x) {
  if(length(x) <= 1) {
    return(x[1]);
  }
  return(HPDinterval(as.mcmc(x),prob = 0.95)[1])
}

summarise_hpd_upper <- function(x) {
  if(length(x) <= 1) {
    return(x[1]);
  }
  return(HPDinterval(as.mcmc(x),prob = 0.95)[2])
}

tmrca_hemi <- rbind(by_tmrca_hemi)
tmrca_hemi1 <- tmrca_hemi %>%
  group_by(virus, season) %>%
  summarise(tmrca_mean = mean(tmrca),
            tmrca_upper = summarise_hpd_upper(tmrca),
            tmrca_lower = summarise_hpd_lower(tmrca))

turnover_by3 <- tmrca_hemi1 %>% filter(virus == "by") %>%
  mutate(time = as.numeric(sapply(str_split(season, "_"), function(x) x[3])))
turnover_by3$y_tmrca <- decimal_date(as.Date("2020-01-01")) - turnover_by3$tmrca_mean                                     
turnover_by3$y_tmrca_lower <- decimal_date(as.Date("2020-01-01")) - turnover_by3$tmrca_lower  
turnover_by3$y_tmrca_upper <- decimal_date(as.Date("2020-01-01")) - turnover_by3$tmrca_upper  
turnover_by3 <- turnover_by3 %>% filter(time >= 2013) %>% filter(time <= 2020)
ggplot(data = turnover_by3) +
  annotate("rect", xmin = 2020.085,xmax = 2021.578,
           ymin = 0,ymax = 18,alpha = 0.2,fill = colors[5])+
  geom_point(aes(x = time, y = y_tmrca), color = colors1[9])+
  geom_errorbar(aes(x = time, ymin = y_tmrca_lower, ymax = y_tmrca_upper),width = 0)+
  scale_y_continuous("Time of MRCA\n(years before Jan 2020)", limits = c(0,18),breaks = seq(0,18,3), expand = c(0,0))+
  scale_x_continuous("Northern and Southern Hemisphere winter season",expand = c(0.01,0),breaks = seq(2013,2023.5,0.5),
                     limits = c(2012.497,2023.575),
                     labels = rev(c("2023","2022/2023","2022","2021/2022","2021","2020/2021","2020",
                                "2019/2020","2019","2018/2019","2018",
                                "2017/2018","2017","2016/2017","2016","2015/2016","2015",
                                "2014/2015","2014","2013/2014","2013","2012/2013")))+
  theme_bw()+
  theme(
    axis.text.x = element_text(angle = 45,hjust =1,vjust = 1),
    # axis.title.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    legend.background = element_blank())+
  annotate("segment", x = 2020, y = 0, xend = 2013, yend = 7, linetype = 2)+
  labs(subtitle = "B/Yamagata", tag = "d")-> p7

#==4. Genetic diversity
diver_h1n1 <- read.delim("//wsl.localhost/Ubuntu-22.04/home/zhiyuanchen/PACT-0.9.5/h1n1_diversity/out.skylines") %>%
  mutate(virus_type = "H1N1pdm09")
diver_h3n2 <- read.delim("//wsl.localhost/Ubuntu-22.04/home/zhiyuanchen/PACT-0.9.5/h3n2_diversity/out.skylines") %>%
  mutate(virus_type = "H3N2")
diver_bv <- read.delim("//wsl.localhost/Ubuntu-22.04/home/zhiyuanchen/PACT-0.9.5/bv_diversity/out.skylines") %>%
  mutate(virus_type = "B/Victoria")
diver_by <- read.delim("//wsl.localhost/Ubuntu-22.04/home/zhiyuanchen/PACT-0.9.5/by_diversity/out.skylines") %>%
  mutate(virus_type = "B/Yamagata")

diver <- rbind(diver_h1n1, diver_h3n2, diver_bv, diver_by) %>% 
  filter(!is.nan(mean)) %>%
  mutate(date = decimal2Date(time))

factor(diver$virus_type, levels = unique(diver$virus_type)) -> diver$virus_type
ggplot(data = diver[diver$date >= as.Date("2011-01-01"),]) +
  annotate("rect", xmin = as.Date("2020-02-01"),xmax = as.Date("2021-07-31"),
           ymin = 0,ymax = 12,alpha = 0.2,fill = colors[5])+
  geom_ribbon(aes(x = date, ymin = lower, ymax = upper, fill = virus_type), alpha = 0.2)+
  geom_line(aes(x = date, y = mean, color = virus_type, group = virus_type))+
  geom_point(aes(x = date, y = mean, color = virus_type))+
  scale_x_date(expand = c(0,0),date_labels = "%Y",date_breaks = "1 year",
               limits = c(as.Date("2012-07-01"),as.Date("2023-07-30")))+
  theme_bw()+
  theme(legend.position = "none",
        legend.background = element_blank())+
  labs(x = "Year", y = "Mean pairwise diversity", tag = "c")+
  scale_y_continuous(limits = c(0,12), expand = c(0,0))+
  scale_color_manual("",values = colors1[c(1:3,9)])+
  scale_fill_manual("",values = colors1[c(1:3,9)]) -> fig_c

#==5.Selection pressure
ratio <- read_xlsx("../phylogeography/4.1post-analyses/nature_selection/ratio.xlsx")
h1n1_slac <- read.csv("../phylogeography/4.1post-analyses/nature_selection/datamonkey-table-h1n1-slac.csv")
h3n2_slac <- read.csv("../phylogeography/4.1post-analyses/nature_selection/datamonkey-table-h3n2-slac.csv")
bv_slac <- read.csv("../phylogeography/4.1post-analyses/nature_selection/datamonkey-table-bv-slac.csv")
by_slac <- read.csv("../phylogeography/4.1post-analyses/nature_selection/datamonkey-table-by-slac.csv")

h1n1_fubar <- read.csv("../phylogeography/4.1post-analyses/nature_selection/datamonkey-table-h1n1-fubar.csv")
h3n2_fubar <- read.csv("../phylogeography/4.1post-analyses/nature_selection/datamonkey-table-h3n2-fubar.csv")
bv_fubar <- read.csv("../phylogeography/4.1post-analyses/nature_selection/datamonkey-table-bv-fubar.csv")
by_fubar <- read.csv("../phylogeography/4.1post-analyses/nature_selection/datamonkey-table-by-fubar.csv")

n1 <- h1n1_slac$Site[h1n1_slac$P..dN.dS...1. < 0.1]
n2 <-h3n2_slac$Site[h3n2_slac$P..dN.dS...1. < 0.1]
n3 <-bv_slac$Site[bv_slac$P..dN.dS...1. < 0.1]
n4 <-by_slac$Site[by_slac$P..dN.dS...1. < 0.1]

m1 <-h1n1_fubar$Site[h1n1_fubar$Prob..alpha...beta...1 > 0.8]
m2 <-h3n2_fubar$Site[h3n2_fubar$Prob..alpha...beta...1 > 0.8]
m3 <-bv_fubar$Site[bv_fubar$Prob..alpha...beta...1 > 0.8]
m4 <-by_fubar$Site[by_fubar$Prob..alpha...beta...1 > 0.8]

ratio$num_pos[1] <- length(n1) + length(m1) - length(unique(c(n1, m1)))
ratio$num_pos[2] <- length(n2) + length(m2) - length(unique(c(n2, m2)))
ratio$num_pos[3] <- length(n3) + length(m3) - length(unique(c(n3, m3)))
ratio$num_pos[4] <- length(n4) + length(m4) - length(unique(c(n4, m4)))

factor(ratio$virus, levels = unique(ratio$virus)) -> ratio$virus
ggplot(data = ratio) +
  geom_bar(aes(x = virus, y = type), stat = "identity", width = 0.3, fill = colors1[c(1:3,9)])+
  geom_line(aes(x = virus, y = num_pos/50, group = 1), linetype = 2)+
  geom_point(aes(x = virus, y = num_pos/50))+
  scale_y_continuous(expression("Selection pressure, d"["N"]/"d"["S"]),
                     sec.axis = sec_axis(trans=~.* 50,
                                         name="No. of sites under\npositive selection"))+
  theme_bw()+
  theme(axis.title.x = element_blank())+
  labs(tag = "e") -> fig_e

#==output
pdf("output/Fig6.pdf", width = 10, height = 10)
p1 -> part1
((p2/fig_c/p7/fig_e)) -> part2
viewport(x = 0, y = 0, width = 0.4, height = 1, just = c("left", "bottom")) -> vp1
viewport(x = 0.4, y = 0, width = 0.6, height = 1, just = c("left", "bottom")) -> vp2
viewport(x = 0.02, y = 0.6, width = 0.19, height = 0.33, just = c("left", "bottom")) -> vp3
print(part1,vp = vp1)
print(part2, vp = vp2)
print(p1_1, vp = vp3)
# (((p2/p1)+ plot_layout(heights = c(0.3,0.7)))|((p3/p4/p5/p6) + plot_layout(guides = "collect")&theme(legend.position = "bottom")))+plot_layout(widths = c(0.4,0.6))+plot_annotation(tag_levels = "a")
dev.off()
