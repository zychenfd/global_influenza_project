setwd("C:/Users/zyche/Nutstore/1/Evolution_study/Oxford Visiting/Influenza project/code")

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

#==define date==
date_match <- data.frame(date = seq(as.Date("1995-01-02"),as.Date("2023-12-25"),7))
date_match$ISO_YEAR <- isoyear(date_match$date)
date_match$ISO_WEEK <- isoweek(date_match$date)

#==h1n1pdm09&h3n2==
#1.read metadata
#GISAID
gisaid_h1n1_meta <- read.csv("../../../flu_data/h1n1/h1n1_ha_meta_high.csv")
gisaid_h1n1_meta <- gisaid_h1n1_meta[,-12]
names(gisaid_h1n1_meta)[c(7,24)] <- c("alpha.3","seqName")
gisaid_h3n2_meta <- read.csv("../../../flu_data/h3n2/h3n2_ha_meta_high.csv")
gisaid_h3n2_meta <- gisaid_h3n2_meta[,-12]
names(gisaid_h3n2_meta)[c(7,24)] <- c("alpha.3","seqName")
#NCBI
ncbi_h1n1_meta <- read.csv("../../../flu_data/h1n1/ncbi_ha_meta_h1n1_high.csv")
ncbi_h3n2_meta <- read.csv("../../../flu_data/h3n2/ncbi_ha_meta_h3n2_high.csv")
#2.read sequence
gisaid_h1n1_seq <- readDNAStringSet("../../../flu_data/h1n1/h1n1.aligned_high_qua.fasta")
gisaid_h3n2_seq <- readDNAStringSet("../../../flu_data/h3n2/h3n2.aligned_high_qua.fasta")
ncbi_h1n1_seq <- readDNAStringSet("../../../flu_data/h1n1/ncbi_seq_h1n1_ha_align_high.fasta")
ncbi_h3n2_seq <- readDNAStringSet("../../../flu_data/h3n2/ncbi_seq_h3n2_ha_align_high.fasta")
#combine metadata and sequence
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

# write.csv(h1n1_meta, "../../../flu_data/Final_data/h1n1_meta.csv", row.names = F)
# write.csv(h3n2_meta, "../../../flu_data/Final_data/h3n2_meta.csv", row.names = F)
# writeXStringSet(h1n1_seq, "../../../flu_data/Final_data/h1n1_seq.fasta", format="fasta")
# writeXStringSet(h3n2_seq, "../../../flu_data/Final_data/h3n2_seq.fasta", format="fasta")

#==B/Victoria&B/Yamagata==
#1.read metadata
#GISAID
gisaid_bv_meta <- read.csv("../../../flu_data/bv_ha_meta_high.csv")
gisaid_bv_meta <- gisaid_bv_meta[,-12]
names(gisaid_bv_meta)[c(7,23)] <- c("alpha.3","seqName")
gisaid_by_meta <- read.csv("../../../flu_data/by_ha_meta_high.csv")
gisaid_by_meta <- gisaid_by_meta[,-12]
names(gisaid_by_meta)[c(7,23)] <- c("alpha.3","seqName")
#NCBI
ncbi_bv_meta <- read.csv("../../../flu_data/ncbi_ha_meta_bv_high.csv") %>% select(names(gisaid_bv_meta))
ncbi_by_meta <- read.csv("../../../flu_data/ncbi_ha_meta_by_high.csv") %>% select(names(gisaid_by_meta))
#2.read sequence
gisaid_bv_seq <- readDNAStringSet("../../../flu_data/bv_ha/nextclade.aligned_high_qua.fasta")
gisaid_by_seq <- readDNAStringSet("../../../flu_data/by_ha/nextclade.aligned_high_qua.fasta")
ncbi_bv_seq <- readDNAStringSet("../../../flu_data/ncbi_seq_bv_ha_align_high.fasta")
ncbi_by_seq <- readDNAStringSet("../../../flu_data/ncbi_seq_by_ha_align_high.fasta")
#combine metadata and sequence
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

# write.csv(bv_meta, "../../../flu_data/Final_data/BV_meta.csv", row.names = F)
# write.csv(by_meta, "../../../flu_data/Final_data/BY_meta.csv", row.names = F)
# writeXStringSet(bv_seq, "../../../flu_data/Final_data/BV_seq.fasta", format="fasta")
# writeXStringSet(by_seq, "../../../flu_data/Final_data/BY_seq.fasta", format="fasta")

#==Epi data from FluNet==
epi <- read.csv("../data/epi_data/VIW_FNT.csv") %>% 
  select("COUNTRY_AREA_TERRITORY","COUNTRY_CODE","ISO2",
         "ISO_YEAR","ISO_WEEK","SPEC_RECEIVED_NB","SPEC_PROCESSED_NB",
         "AH1N12009","AH3","ANOTSUBTYPED","ANOTSUBTYPABLE","INF_A",
         "BVIC_NODEL","BYAM","BNOTDETERMINED","INF_B","INF_NEGATIVE") %>%
  filter(ISO_YEAR >= 2010) %>% # 119314
  filter(!(is.na(SPEC_PROCESSED_NB) & is.na(INF_NEGATIVE))) %>%
  filter(!(SPEC_PROCESSED_NB == 0 & is.na(INF_NEGATIVE))) %>%
  filter(!(SPEC_PROCESSED_NB == 0 & INF_NEGATIVE == 0)) %>%
  filter(!(is.na(INF_A) & is.na(INF_B) & is.na(INF_NEGATIVE)))
epi$COUNTRY_CODE[str_detect(epi$COUNTRY_AREA_TERRITORY,"United Kingdom")] <- "GBR"

epi$INF_NEGATIVE[is.na(epi$INF_NEGATIVE)] <- 
  epi$SPEC_PROCESSED_NB[is.na(epi$INF_NEGATIVE)]-
  epi$INF_A[is.na(epi$INF_NEGATIVE)]-
  epi$INF_B[is.na(epi$INF_NEGATIVE)]
epi <- epi %>% filter(INF_NEGATIVE >= -5)

epi$AH1N12009[is.na(epi$AH1N12009)] <- 0
epi$AH3[is.na(epi$AH3)] <- 0
epi$ANOTSUBTYPED[is.na(epi$ANOTSUBTYPED)] <- 0
epi$ANOTSUBTYPABLE[is.na(epi$ANOTSUBTYPABLE)] <- 0
epi$INF_A[is.na(epi$INF_A)] <- 0
epi$BVIC_NODEL[is.na(epi$BVIC_NODEL)] <- 0
epi$BYAM[is.na(epi$BYAM)] <- 0
epi$BNOTDETERMINED[is.na(epi$BNOTDETERMINED)] <- 0
epi$INF_B[is.na(epi$INF_B)] <- 0
epi$INF_NEGATIVE[is.na(epi$INF_NEGATIVE)] <- 0

epi$SPEC_PROCESSED_NB[is.na(epi$SPEC_PROCESSED_NB) | epi$SPEC_PROCESSED_NB == 0] <- 
  epi$INF_A[is.na(epi$SPEC_PROCESSED_NB) | epi$SPEC_PROCESSED_NB == 0] +
  epi$INF_B[is.na(epi$SPEC_PROCESSED_NB) | epi$SPEC_PROCESSED_NB == 0] +
  epi$INF_NEGATIVE[is.na(epi$SPEC_PROCESSED_NB) | epi$SPEC_PROCESSED_NB == 0]

#logical check
table(epi$INF_A >= c(epi$AH1N12009 + epi$AH3 + epi$ANOTSUBTYPED + epi$ANOTSUBTYPABLE)) # T
table(epi$INF_B == c(epi$BVIC_NODEL + epi$BYAM + epi$BNOTDETERMINED)) # T
# table(epi$SPEC_RECEIVED_NB >= epi$SPEC_PROCESSED_NB) #logical error 1
table(epi$SPEC_PROCESSED_NB >= c(epi$INF_A+epi$INF_B)) #infected with fluA and fluB at the same time
table(epi$SPEC_PROCESSED_NB == c(epi$INF_A+epi$INF_B+epi$INF_NEGATIVE))

epi <- epi %>%
  filter(SPEC_PROCESSED_NB - INF_A - INF_B - INF_NEGATIVE >= -5) %>%
  filter(SPEC_PROCESSED_NB - INF_A - INF_B - INF_NEGATIVE <= 0)

epi <- left_join(epi, date_match)
epi$A_notype <- epi$ANOTSUBTYPED + epi$ANOTSUBTYPABLE
epi$A_other <- epi$INF_A - epi$AH1N12009 - epi$AH3 - epi$A_notype
epi1 <- epi %>% 
  filter(A_other >= 0) %>%
  group_by(COUNTRY_CODE,
           ISO_YEAR,ISO_WEEK,date) %>%
  summarise(SPEC_PROCESSED_NB = sum(SPEC_PROCESSED_NB),
            AH1N12009 = sum(AH1N12009),
            AH3 = sum(AH3),
            A_other = sum(A_other),
            A_notype = sum(A_notype),
            INF_A = sum(INF_A),
            BVIC_NODEL = sum(BVIC_NODEL),
            BYAM = sum(BYAM),
            BNOTDETERMINED =  sum(BNOTDETERMINED),
            INF_B = sum(INF_B),
            INF_NEGATIVE = sum(INF_NEGATIVE)) %>%
  filter(date >= as.Date("2010-01-04"))

len <- length(seq(as.Date("2010-01-04"),as.Date("2023-08-28"),7))
epi2 <- data.frame(COUNTRY_CODE = rep(unique(epi$COUNTRY_CODE),each = len),
                   date = rep(seq(as.Date("2010-01-04"),as.Date("2023-08-28"),7),
                              length = length(unique(epi$COUNTRY_CODE))*len ))
epi2 <- left_join(epi2, epi1)

epi2$h1n1_prop <- NA
epi2$h1n1_prop_LL <- NA
epi2$h1n1_prop_UL <- NA
epi2$h3n2_prop <- NA
epi2$h3n2_prop_LL <- NA
epi2$h3n2_prop_UL <- NA
for (i in 1:nrow(epi2)) {
  # i = 1
  if (is.na(epi2$A_notype[i])) {
    epi2$h1n1_prop[i] <- NA
    epi2$h3n2_prop[i] <- NA
  } else {
    if (epi2$A_notype[i] >= 1) {
      tmp <- epi2 %>% filter(COUNTRY_CODE == epi2$COUNTRY_CODE[i] &
                               date <= (epi2$date[i] + 0) & 
                               date >= (epi2$date[i] - 0))
      tmp1 <- tmp %>% summarise(h1n1 = sum(AH1N12009, na.rm = T),
                                h3n2 = sum(AH3,  na.rm = T),
                                other = sum(A_other,  na.rm = T)) %>%
        mutate(h1n1_prop = (h1n1+1)/(h1n1+h3n2+other+2),
               h1n1_prop_LL = qbeta(0.025, 1+h1n1, h3n2+other + 1),
               h1n1_prop_UL = qbeta(0.975, 1+h1n1, h3n2+other + 1),
               h3n2_prop = (h3n2+1)/(h1n1+h3n2+other+2),
               h3n2_prop_LL = qbeta(0.025, 1+h3n2, h1n1+other + 1),
               h3n2_prop_UL = qbeta(0.975, 1+h3n2, h1n1+other + 1))
      epi2$h1n1_prop[i] <- tmp1$h1n1_prop
      epi2$h3n2_prop[i] <- tmp1$h3n2_prop
      epi2$h1n1_prop_LL[i] <- tmp1$h1n1_prop_LL
      epi2$h1n1_prop_UL[i] <- tmp1$h1n1_prop_UL
      epi2$h3n2_prop_LL[i] <- tmp1$h3n2_prop_LL
      epi2$h3n2_prop_UL[i] <- tmp1$h3n2_prop_UL
    } else {
      epi2$h1n1_prop[i] <- NA
      epi2$h3n2_prop[i] <- NA
    }
  }
}

epi2$bv_prop <- NA
epi2$bv_prop_LL <- NA
epi2$bv_prop_UL <- NA
for (i in 1:nrow(epi2)) {
  # i = 2
  if (is.na(epi2$BNOTDETERMINED[i])) {
    epi2$bv_prop[i] <- NA
  } else {
    if (epi2$BNOTDETERMINED[i] >= 1) {
      tmp <- epi2 %>% filter(COUNTRY_CODE == epi2$COUNTRY_CODE[i] &
                               date <= (epi2$date[i] + 0) & 
                               date >= (epi2$date[i] - 0))
      tmp1 <- tmp %>% summarise(bv = sum(BVIC_NODEL, na.rm = T),
                                by = sum(BYAM,  na.rm = T)) %>%
        mutate(bv_prop = (bv+1)/(bv+by+2),
               bv_prop_LL = qbeta(0.025, 1+bv, by + 1),
               bv_prop_UL = qbeta(0.975, 1+bv, by + 1),)
      epi2$bv_prop[i] <- tmp1$bv_prop
      epi2$bv_prop_LL[i] <- tmp1$bv_prop_LL
      epi2$bv_prop_UL[i] <- tmp1$bv_prop_UL
    } else {
      epi2$bv_prop[i] <- NA
    }
  }
}

epi2$by_prop <- NA
epi2$by_prop_LL <- NA
epi2$by_prop_UL <- NA
for (i in 1:nrow(epi2)) {
  # i = 2
  if (is.na(epi2$BNOTDETERMINED[i])) {
    epi2$by_prop[i] <- NA
  } else {
    if (epi2$BNOTDETERMINED[i] >= 1) {
      tmp <- epi2 %>% filter(COUNTRY_CODE == epi2$COUNTRY_CODE[i] &
                               date <= (epi2$date[i] + 0) & 
                               date >= (epi2$date[i] - 0))
      tmp1 <- tmp %>% summarise(bv = sum(BVIC_NODEL, na.rm = T),
                                by = sum(BYAM,  na.rm = T)) %>%
        mutate(by_prop = (by+1)/(bv+by+2),
               by_prop_LL = qbeta(0.025, 1+by, bv + 1),
               by_prop_UL = qbeta(0.975, 1+by, bv + 1),)
      epi2$by_prop[i] <- tmp1$by_prop
      epi2$by_prop_LL[i] <- tmp1$by_prop_LL
      epi2$by_prop_UL[i] <- tmp1$by_prop_UL
    } else {
      epi2$by_prop[i] <- NA
    }
  }
}

nrow(epi2 %>%
  filter(A_notype > 0) %>%
  filter(is.na(h1n1_prop)))

#==
epi_final <- epi2
epi_final <- epi_final %>% filter(!is.na(SPEC_PROCESSED_NB))
epi_final$h1n1_no <- epi_final$AH1N12009 + epi_final$A_notype * epi_final$h1n1_prop
epi_final$h1n1_no_LL <- epi_final$AH1N12009 + epi_final$A_notype * epi_final$h1n1_prop_LL
epi_final$h1n1_no_UL <- epi_final$AH1N12009 + epi_final$A_notype * epi_final$h1n1_prop_UL
epi_final$h1n1_no[is.na(epi_final$h1n1_no)] <- epi_final$AH1N12009[is.na(epi_final$h1n1_no)]
epi_final$h1n1_no_LL[is.na(epi_final$h1n1_no_LL)] <- epi_final$AH1N12009[is.na(epi_final$h1n1_no_LL)]
epi_final$h1n1_no_UL[is.na(epi_final$h1n1_no_UL)] <- epi_final$AH1N12009[is.na(epi_final$h1n1_no_UL)]

epi_final$h3n2_no <- epi_final$AH3 + epi_final$A_notype * epi_final$h3n2_prop
epi_final$h3n2_no_LL <- epi_final$AH3 + epi_final$A_notype * epi_final$h3n2_prop_LL
epi_final$h3n2_no_UL <- epi_final$AH3 + epi_final$A_notype * epi_final$h3n2_prop_UL
epi_final$h3n2_no[is.na(epi_final$h3n2_no)] <- epi_final$AH3[is.na(epi_final$h3n2_no)]
epi_final$h3n2_no_LL[is.na(epi_final$h3n2_no_LL)] <- epi_final$AH3[is.na(epi_final$h3n2_no_LL)]
epi_final$h3n2_no_UL[is.na(epi_final$h3n2_no_UL)] <- epi_final$AH3[is.na(epi_final$h3n2_no_UL)]

epi_final$bv_no <- epi_final$BVIC_NODEL + epi_final$BNOTDETERMINED * epi_final$bv_prop
epi_final$bv_no_LL <- epi_final$BVIC_NODEL + epi_final$BNOTDETERMINED * epi_final$bv_prop_LL
epi_final$bv_no_UL <- epi_final$BVIC_NODEL + epi_final$BNOTDETERMINED * epi_final$bv_prop_UL
epi_final$bv_no[is.na(epi_final$bv_no)] <- epi_final$BVIC_NODEL[is.na(epi_final$bv_no)]
epi_final$bv_no_LL[is.na(epi_final$bv_no_LL)] <- epi_final$BVIC_NODEL[is.na(epi_final$bv_no_LL)]
epi_final$bv_no_UL[is.na(epi_final$bv_no_UL)] <- epi_final$BVIC_NODEL[is.na(epi_final$bv_no_UL)]

epi_final$by_no <- epi_final$BYAM + epi_final$BNOTDETERMINED * epi_final$by_prop
epi_final$by_no_LL <- epi_final$BYAM + epi_final$BNOTDETERMINED * epi_final$by_prop_LL
epi_final$by_no_UL <- epi_final$BYAM + epi_final$BNOTDETERMINED * epi_final$by_prop_UL
epi_final$by_no[is.na(epi_final$by_no)] <- epi_final$BYAM[is.na(epi_final$by_no)]
epi_final$by_no_LL[is.na(epi_final$by_no_LL)] <- epi_final$BYAM[is.na(epi_final$by_no_LL)]
epi_final$by_no_UL[is.na(epi_final$by_no_UL)] <- epi_final$BYAM[is.na(epi_final$by_no_UL)]

epi_final$diff <- epi_final$SPEC_PROCESSED_NB - epi_final$INF_A - epi_final$INF_B- epi_final$INF_NEGATIVE

#add region information
world_region <- read.csv("C:/Users/zyche/Nutstore/1/Evolution_study/COVID/China/data/list_country.csv")
who_country <- read_xlsx("../data/WHO_country.xlsx")
epi_final <- left_join(epi_final, world_region[,c(3,6:8)],by = c("COUNTRY_CODE"="alpha.3"))
epi_final$region_final <- epi_final$region
epi_final$region_final[epi_final$region_final == "Europe" & epi_final$COUNTRY_CODE == "RUS"] <- "Russia"
epi_final$region_final[epi_final$COUNTRY_CODE == "XKX"] <- "Europe"
epi_final$region_final[epi_final$region_final == "Americas"] <- "Northern America"
epi_final$region_final[epi_final$intermediate.region == "South America"] <- "Southern America"
epi_final$region_final[epi_final$region_final == "Asia"] <- epi_final$sub.region[epi_final$region_final == "Asia"]
epi_final$region_final[epi_final$COUNTRY_CODE %in% c("CHN","HKG","TWN","MAC")] <- "China"
epi_final$region_final[epi_final$COUNTRY_CODE %in% c("KOR","JPN")] <- "Japan/Korea"
table(epi_final$region_final)
unique(epi_final$COUNTRY_CODE)[which(unique(epi_final$COUNTRY_CODE) %in% c(who_country$ISO3,"HKG","TWN","MAC") == F)]
epi_final <- epi_final %>%
  filter(epi_final$COUNTRY_CODE %in% c(who_country$ISO3,"HKG","TWN","MAC"))

#by region
epi_glo_flu_region <- epi_final %>% 
  group_by(ISO_YEAR, ISO_WEEK, region_final) %>% 
  summarise(test = sum(SPEC_PROCESSED_NB),
            h1n1_num = sum(h1n1_no,na.rm = T),
            h1n1_num_LL = sum(h1n1_no_LL,na.rm = T),
            h1n1_num_UL = sum(h1n1_no_UL,na.rm = T),
            h3n2_num = sum(h3n2_no,na.rm = T),
            h3n2_num_LL = sum(h3n2_no_LL,na.rm = T),
            h3n2_num_UL = sum(h3n2_no_UL,na.rm = T),
            BV_num = sum(bv_no, na.rm = T),
            BV_num_LL = sum(bv_no_LL, na.rm = T),
            BV_num_UL = sum(bv_no_UL, na.rm = T),
            BY_num = sum(by_no, na.rm = T),
            BY_num_LL = sum(by_no_LL, na.rm = T),
            BY_num_UL = sum(by_no_UL, na.rm = T)) %>%
  # mutate(h1n1_pos = h1n1_num/test,
  #        h1n1_pos_LL = h1n1_num_LL/test,
  #        h1n1_pos_UL = h1n1_num_UL/test,
  #        h3n2_pos = h3n2_num/test,
  #        h3n2_pos_LL = h3n2_num_LL/test,
  #        h3n2_pos_UL = h3n2_num_UL/test,
  #        BV_pos = BV_num/test,
  #        BV_pos_LL = BV_num_LL/test,
  #        BV_pos_UL = BV_num_UL/test,
  #        BY_pos = BY_num/test,
  #        BY_pos_LL = BY_num_LL/test,
  #        BY_pos_UL = BY_num_UL/test) %>%
  left_join(date_match)

#Flu A
epi_glo_flua <- epi_final %>% 
  group_by(ISO_YEAR, ISO_WEEK) %>% 
  summarise(test = sum(SPEC_PROCESSED_NB),
            h1n1_num = sum(h1n1_no,na.rm = T),
            h1n1_num_LL = sum(h1n1_no_LL,na.rm = T),
            h1n1_num_UL = sum(h1n1_no_UL,na.rm = T),
            h3n2_num = sum(h3n2_no,na.rm = T),
            h3n2_num_LL = sum(h3n2_no_LL,na.rm = T),
            h3n2_num_UL = sum(h3n2_no_UL,na.rm = T)) %>%
  # mutate(h1n1_pos = h1n1_num/test,
  #        h1n1_pos_LL = h1n1_num_LL/test,
  #        h1n1_pos_UL = h1n1_num_UL/test,
  #        h3n2_pos = h3n2_num/test,
  #        h3n2_pos_LL = h3n2_num_LL/test,
  #        h3n2_pos_UL = h3n2_num_UL/test) %>%
  left_join(date_match)

#Flu B
epi_glo_flub <- epi_final %>% 
  group_by(ISO_YEAR, ISO_WEEK) %>% 
  summarise(test = sum(SPEC_PROCESSED_NB),
            BV_num = sum(bv_no, na.rm = T),
            BV_num_LL = sum(bv_no_LL, na.rm = T),
            BV_num_UL = sum(bv_no_UL, na.rm = T),
            BY_num = sum(by_no, na.rm = T),
            BY_num_LL = sum(by_no_LL, na.rm = T),
            BY_num_UL = sum(by_no_UL, na.rm = T)) %>%
  # mutate(BV_pos = BV_num/test,
  #        BV_pos_LL = BV_num_LL/test,
  #        BV_pos_UL = BV_num_UL/test,
  #        BY_pos = BY_num/test,
  #        BY_pos_LL = BY_num_LL/test,
  #        BY_pos_UL = BY_num_UL/test) %>%
  left_join(date_match)

h1n1_meta$date <- as.Date(h1n1_meta$date)
h3n2_meta$date <- as.Date(h3n2_meta$date)
bv_meta$date <- as.Date(bv_meta$date)
by_meta$date <- as.Date(by_meta$date)

#rolling rate
epi_glo_flu_region <- epi_glo_flu_region[order(epi_glo_flu_region$region_final,epi_glo_flu_region$date),]
epi_glo_flu_region1 <- epi_glo_flu_region[,c(3:17)] %>%
  mutate(h1n1_num_roll = rollmean(h1n1_num, k=5, fill=NA, align='center'),
         h1n1_num_roll_LL = rollmean(h1n1_num_LL, k=5, fill=NA, align='center'),
         h1n1_num_roll_UL = rollmean(h1n1_num_UL, k=5, fill=NA, align='center'),
         h3n2_num_roll = rollmean(h3n2_num, k=5, fill=NA, align='center'),
         h3n2_num_roll_LL = rollmean(h3n2_num_LL, k=5, fill=NA, align='center'),
         h3n2_num_roll_UL = rollmean(h3n2_num_UL, k=5, fill=NA, align='center'),
         BV_num_roll = rollmean(BV_num, k=5, fill=NA, align='center'),
         BV_num_roll_LL = rollmean(BV_num_LL, k=5, fill=NA, align='center'),
         BV_num_roll_UL = rollmean(BV_num_UL, k=5, fill=NA, align='center'),
         BY_num_roll = rollmean(BY_num, k=5, fill=NA, align='center'),
         BY_num_roll_LL = rollmean(BY_num_LL, k=5, fill=NA, align='center'),
         BY_num_roll_UL = rollmean(BY_num_UL, k=5, fill=NA, align='center'),
         test_roll = rollmean(test, k=5, fill=NA, align='center'))
epi_glo_flu_region1$diff <- c(0,diff.Date(epi_glo_flu_region1$date))

epi_glo_flua1 <- epi_glo_flua[,c(3:10)] %>%
  mutate(h1n1_num_roll = rollmean(h1n1_num, k=5, fill=NA, align='center'),
         h1n1_num_roll_LL = rollmean(h1n1_num_LL, k=5, fill=NA, align='center'),
         h1n1_num_roll_UL = rollmean(h1n1_num_UL, k=5, fill=NA, align='center'),
         h3n2_num_roll = rollmean(h3n2_num, k=5, fill=NA, align='center'),
         h3n2_num_roll_LL = rollmean(h3n2_num_LL, k=5, fill=NA, align='center'),
         h3n2_num_roll_UL = rollmean(h3n2_num_UL, k=5, fill=NA, align='center'),
         test_roll = rollmean(test, k=5, fill=NA, align='center'))
epi_glo_flua1$diff <- c(0,diff.Date(epi_glo_flua1$date))

epi_glo_flub1 <- epi_glo_flub[,c(3:10)] %>%
  mutate(BV_num_roll = rollmean(BV_num, k=5, fill=NA, align='center'),
         BV_num_roll_LL = rollmean(BV_num_LL, k=5, fill=NA, align='center'),
         BV_num_roll_UL = rollmean(BV_num_UL, k=5, fill=NA, align='center'),
         BY_num_roll = rollmean(BY_num, k=5, fill=NA, align='center'),
         BY_num_roll_LL = rollmean(BY_num_LL, k=5, fill=NA, align='center'),
         BY_num_roll_UL = rollmean(BY_num_UL, k=5, fill=NA, align='center'),
         test_roll = rollmean(test, k=5, fill=NA, align='center'))
epi_glo_flub1$diff <- c(0,diff.Date(epi_glo_flub1$date))

#save data
# write.csv(epi_final,"../data/epi_data/epi_final_no_roll.csv", row.names = F)
# write.csv(epi_glo_flu_region1,"../data/epi_data/epi_glo_flu_region1_no_roll.csv", row.names = F)
# write.csv(epi_glo_flua1,"../data/epi_data/epi_glo_flua1_no_roll.csv", row.names = F)
# write.csv(epi_glo_flub1,"../data/epi_data/epi_glo_flub1_no_roll.csv", row.names = F)
epi_final <- read.csv("../data/epi_data/epi_final_no_roll.csv") %>% mutate(date = as.Date(date))
epi_glo_flu_region1 <- read.csv("../data/epi_data/epi_glo_flu_region1_no_roll.csv") %>% mutate(date = as.Date(date))
epi_glo_flua1 <- read.csv("../data/epi_data/epi_glo_flua1_no_roll.csv") %>% mutate(date = as.Date(date))
epi_glo_flub1 <- read.csv("../data/epi_data/epi_glo_flub1_no_roll.csv") %>% mutate(date = as.Date(date))

# epi_glo_flua1_stat <- epi_glo_flua1 %>%
#   mutate(rate_h1n1_raw = h1n1_num / test,
#          rate_h3n2_raw = h3n2_num / test,
#          rate_h1n1_roll = h1n1_num_roll / test_roll,
#          rate_h3n2_roll = h3n2_num_roll / test_roll) %>%
#   filter(date >= as.Date("2020-05-01") & date <= as.Date("2021-07-31"))
# epi_glo_flub1_stat <- epi_glo_flub1 %>%
#   mutate(rate_BV_raw = BV_num / test,
#          rate_BY_raw = BY_num / test,
#          rate_BV_roll = BV_num_roll / test_roll,
#          rate_BY_roll = BY_num_roll / test_roll) %>%
#   filter(date >= as.Date("2020-05-01") & date <= as.Date("2021-07-31"))

tmp <- epi_glo_flu_region1 %>% mutate(epi_case = h1n1_num + h3n2_num + BV_num + BY_num) %>%
  group_by(date) %>%
  summarise(epi_case1 = sum(epi_case))
ggplot(data = tmp) + geom_point(aes(x = date, y = epi_case1))+
  scale_x_date(date_breaks = "1 year")

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

#==Postivity rate by region
ggplot(data = epi_glo_flu_region1[!epi_glo_flu_region1$region_final %in% c("Central Asia","Eastern Asia"),]) +
  geom_line(aes(x = date, y = h1n1_num_roll/test_roll, color = region_final))+
  scale_x_date("", date_breaks = "1 year",date_labels = "%Y",
               limits = c(as.Date("2010-01-01"),as.Date("2023-08-31")))+
  scale_y_continuous("H1N1pdm positivity rates") +
  theme_bw() -> x1

ggplot(data = epi_glo_flu_region1[!epi_glo_flu_region1$region_final %in% c("Central Asia","Eastern Asia"),]) +
  geom_line(aes(x = date, y = h3n2_num_roll/test_roll, color = region_final))+
  scale_x_date("", date_breaks = "1 year",date_labels = "%Y",
               limits = c(as.Date("2010-01-01"),as.Date("2023-08-31")))+
  scale_y_continuous("H3N2 positivity rates") +
  theme_bw() -> x2

ggplot(data = epi_glo_flu_region1[!epi_glo_flu_region1$region_final %in% c("Central Asia","Eastern Asia"),]) +
  geom_line(aes(x = date, y = BV_num_roll/test_roll, color = region_final))+
  scale_x_date("", date_breaks = "1 year",date_labels = "%Y",
               limits = c(as.Date("2010-01-01"),as.Date("2023-08-31")))+
  scale_y_continuous("BV positivity rates") +
  theme_bw() -> x3

ggplot(data = epi_glo_flu_region1[!epi_glo_flu_region1$region_final %in% c("Central Asia","Eastern Asia"),]) +
  geom_line(aes(x = date, y = BY_num_roll/test_roll, color = region_final))+
  scale_x_date("", date_breaks = "1 year",date_labels = "%Y",
               limits = c(as.Date("2010-01-01"),as.Date("2023-08-31")))+
  scale_y_continuous("BY positivity rates") +
  theme_bw() -> x4

# tiff(filename = "../output/Postive_rate_by_region.tif",height = 10, width = 10,
#      units = "in", res = 326, compression = "lzw")
# x1/x2/x3/x4+plot_annotation(tag_levels = "a")
# dev.off()

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

###by-regions==
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

library(sf)
library(rgdal)
#global map
world_region <- read.csv("C:/Users/zyche/Nutstore/1/Evolution_study/Oxford Visiting/Influenza project/data/world_region.csv")
nine <- st_read("C:/Users/zyche/Nutstore/1/Evolution_study/Flu_phylogeography/analyses/scripts/data/Geographic/nine.shp")
nine <- st_transform(nine, CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"))
worldmap0 <- st_read("C:/Users/zyche/Nutstore/1/Evolution_study/Flu_phylogeography/analyses/scripts/data/Geographic/World_map.shp")
worldmap <- left_join(worldmap0, world_region[,c(3,4)], by = c("iso_a3" = "alpha.3")) %>% filter(!is.na(region_final))
worldmap$region_final[worldmap$region_final %in% c("China","South China")] <- "China"
worldmap1 <- worldmap[,c(16,17)] %>% group_by(region_final) %>% summarise()
worldmap1 <- st_transform(worldmap1, CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"))
worldmap0 <- st_transform(worldmap0, CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"))

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

library(patchwork)
library(grid)
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


#==transmission level (ratio of positive influenza cases to the total number of specimens processed in each time period)==
epi_glo_flua2 <- epi_glo_flua1 %>%
  filter(date >= as.Date("2018-02-01") & date <= as.Date("2023-07-31")) %>%
  mutate(period = ifelse(date >= as.Date("2021-08-01"),"Post-pandeic","Pandemic")) %>%
  mutate(period = ifelse(date <= as.Date("2020-01-31"),"Pre-pandeic",period)) %>%
  group_by(period) %>%
  summarise(test_total = sum(test))

epi_glo_flua3 <- epi_glo_flua1 %>%
  filter(date >= as.Date("2018-02-01") & date <= as.Date("2023-07-31")) %>%
  mutate(period = ifelse(date >= as.Date("2021-08-01"),"Post-pandeic","Pandemic")) %>%
  mutate(period = ifelse(date <= as.Date("2020-01-31"),"Pre-pandeic",period)) %>%
  left_join(epi_glo_flua2)

ggplot(epi_glo_flua3) +
  annotate("rect", xmin = as.Date("2020-02-01"),xmax = as.Date("2021-07-31"),
           ymin = 0,ymax = 0.0045,alpha = 0.2,fill = colors[5])+
  geom_ribbon(aes(x = date, ymin = h1n1_num_roll_LL/test_total,ymax = h1n1_num_roll_UL/test_total),
              stat = "identity",fill = "lightblue", alpha = 0.8) +
  geom_line(aes(x = date, y = h1n1_num_roll/test_total))+
  scale_x_date("Date", date_breaks = "1 year",date_labels = "%b %Y", limits = c(as.Date("2018-01-01"),as.Date("2023-08-01") ))+
  scale_y_continuous("Transmission level",limits = c(-0.0002, 0.0045), expand = c(0,0))+
  theme_bw()+
  theme(panel.grid.major  = element_blank(),
        axis.title.x = element_blank(),
        plot.margin =  margin(0, 0.5, 0, 0, "cm"),
        panel.grid.minor  = element_blank())+
  labs(subtitle = "H1N1pdm09") -> s1

ggplot(epi_glo_flua3) +
  annotate("rect", xmin = as.Date("2020-02-01"),xmax = as.Date("2021-07-31"),
           ymin = 0,ymax = 0.0045,alpha = 0.2,fill = colors[5])+
  geom_ribbon(aes(x = date, ymin = h3n2_num_roll_LL/test_total,ymax = h3n2_num_roll_UL/test_total),
              stat = "identity",fill = "lightblue", alpha = 0.8) +
  geom_line(aes(x = date, y = h3n2_num_roll/test_total))+
  scale_x_date("Date", date_breaks = "1 year",date_labels = "%b %Y", limits = c(as.Date("2018-01-01"),as.Date("2023-08-01") ))+
  scale_y_continuous("Transmission level",limits = c(-0.0002, 0.0045), expand = c(0,0))+
  theme_bw()+
  theme(panel.grid.major  = element_blank(),
        axis.title.x = element_blank(),
        plot.margin =  margin(0, 0.5, 0, 0, "cm"),
        panel.grid.minor  = element_blank())+
  labs(subtitle = "H3N2") -> s2


epi_glo_flub2 <- epi_glo_flub1 %>%
  filter(date >= as.Date("2018-02-01") & date <= as.Date("2023-07-31")) %>%
  mutate(period = ifelse(date >= as.Date("2021-08-01"),"Post-pandeic","Pandemic")) %>%
  mutate(period = ifelse(date <= as.Date("2020-01-31"),"Pre-pandeic",period)) %>%
  group_by(period) %>%
  summarise(test_total = sum(test))

epi_glo_flub3 <- epi_glo_flub1 %>%
  filter(date >= as.Date("2018-02-01") & date <= as.Date("2023-07-31")) %>%
  mutate(period = ifelse(date >= as.Date("2021-08-01"),"Post-pandeic","Pandemic")) %>%
  mutate(period = ifelse(date <= as.Date("2020-01-31"),"Pre-pandeic",period)) %>%
  left_join(epi_glo_flub2)

ggplot(epi_glo_flub3) +
  annotate("rect", xmin = as.Date("2020-02-01"),xmax = as.Date("2021-07-31"),
           ymin = 0,ymax = 0.0045,alpha = 0.2,fill = colors[5])+
  geom_ribbon(aes(x = date, ymin = BV_num_roll_LL/test_total,ymax = BV_num_roll_UL/test_total),
              stat = "identity",fill = "lightblue", alpha = 0.8) +
  geom_line(aes(x = date, y = BV_num_roll/test_total))+
  scale_x_date("Date", date_breaks = "1 year",date_labels = "%b %Y", limits = c(as.Date("2018-01-01"),as.Date("2023-08-01") ))+
  scale_y_continuous("Transmission level",limits = c(-0.0002, 0.0045), expand = c(0,0))+
  theme_bw()+
  theme(panel.grid.major  = element_blank(),
        axis.title.x = element_blank(),
        plot.margin =  margin(0, 0.5, 0, 0, "cm"),
        panel.grid.minor  = element_blank())+
  labs(subtitle = "B/Victoria") -> s3

ggplot(epi_glo_flub3[epi_glo_flub3$date <= as.Date("2020-03-31"),]) +
  annotate("rect", xmin = as.Date("2020-02-01"),xmax = as.Date("2021-07-31"),
           ymin = 0,ymax = 0.0045,alpha = 0.2,fill = colors[5])+
  geom_ribbon(aes(x = date, ymin = BY_num_roll_LL/test_total,ymax = BY_num_roll_UL/test_total),
              stat = "identity",fill = "lightblue", alpha = 0.8) +
  geom_line(aes(x = date, y = BY_num_roll/test_total))+
  scale_x_date("Date", date_breaks = "1 year",date_labels = "%b %Y", limits = c(as.Date("2018-01-01"),as.Date("2023-08-01") ))+
  scale_y_continuous("Transmission level",limits = c(-0.0002, 0.0045), expand = c(0,0))+
  theme_bw()+
  theme(panel.grid.major  = element_blank(),
        plot.margin =  margin(0, 0.5, 0, 0, "cm"),
        panel.grid.minor  = element_blank())+
  labs(subtitle = "B/Yamagata") -> s4

pdf("../output/FigS_activity_level.pdf",height = 7, width = 7)
s1/s2/s3/s4+plot_annotation(tag_levels = "a")
dev.off()

svg("../output/FigS_activity_level.svg",height = 7, width = 7)
s1/s2/s3+plot_annotation(tag_levels = "a")
dev.off()
