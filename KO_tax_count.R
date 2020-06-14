#R v.3.6.1
setwd("/Users/Winni/Desktop/Mueller_Lab/2019_Ischia/Metagenomes/KEGG")
library(tidyr);library(dplyr)

KO_family <- read.table(file="KO_family_R.txt", header=TRUE, sep="\t")
head(KO_family)

#removing empty rows
KO_family2 <- KO_family[-which(KO_family$Family ==""), ]
KO_family3 <- KO_family2[-which(KO_family2$Family=="Unclassified"), ]
head(KO_family3)

#counting the number of occurrences for family 
KO_family_grouped <- KO_family3 %>%
  group_by(Pathway, Treatment, Family) %>%
  summarise(n=n())
head(KO_family_grouped)

write.table(KO_family_grouped2, "KO_family_grouped.txt", quote=FALSE, sep="\t")
head(KO_family_grouped)

#getting rel abund info
KO_family_grouped2 <- KO_family_grouped %>%
  group_by(Pathway, Treatment) %>%
  mutate(relAbundByPath_Treat = (n/ sum(n))*100)

setwd("/Users/Winni/Desktop/Mueller_Lab/2019_Ischia/Metagenomes/")
KO_genus <-read.table(file="KO_grouped_genus.txt", header=TRUE, sep="\t")
KO_genus2 <- KO_genus[-which(KO_genus$Genus ==""), ]
KO_genus3 <- KO_genus2[-which(KO_genus2$Genus=="Unclassified"), ]

KO_genus_grouped <- KO_genus3 %>%
  group_by(Pathway, Treatment, Genus) %>%
  summarise(n=n())

KO_genus_grouped2 <- KO_genus_grouped %>%
  group_by(Pathway, Treatment) %>%
  mutate(relAbundByPath_Treat = (n/sum(n))*100)

write.table(KO_genus_grouped2, "KO_genus_grouped.txt", quote=FALSE, sep="\t")
