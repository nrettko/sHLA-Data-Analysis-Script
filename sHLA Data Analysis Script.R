## Purose: Filter sHLA Mass Spec Data for HLA-Associated Peptides via NetMHC 4.0 Server - DTU
## Author: Nick Rettko (Nicholas.Rettko@ucsf.edu)
## Updated: February 15, 2021
##Code will only work if you save Search Compare results from Protein Prospector as a .txt file and later convert to .csv

install.packages("data.table")
install.packages("stringi")
install.packages("stringr")
install.packages("seqinr")
install.packages("ggplot2")

library(data.table)
library(stringi)
library(stringr)
library(seqinr)

#Set working directory. This is the folder where everything will take place and it will look for other folders.
setwd("~/Documents/Wells Lab/sHLA Results/MCF7/")

Rawdata <- read.csv(file="~/Documents/Wells Lab/sHLA Results/MCF7/210724_sHLA_MCF7_Senescence_Growing.db.MCF7 Senescence 1.peptides.csv", stringsAsFactors = FALSE)

#Removes everything but the UniProt ID, Peptide Sequence, and Protein Name and counts the Peptide Length
RawTable <- Rawdata[,-c(2,5:15,18:19)]
#RawTable <- Rawdata[-c(1),-c(2,5:16,19:20)]
UniProtID <- gsub("\\|.*","",as.character(unlist(c(RawTable[,c(4)]))))
Peptide_Sequence <- as.character(unlist(c(RawTable[,c(1)])))
#Peptide_Sequence <- str_replace_all(Peptide_Sequence, "[^[:alnum:]]", "")
#Peptide_Sequence <- str_replace_all(Peptide_Sequence, "[[:digit:]]+", "")
Protein_Name <- gsub('^.','',sub("^[^|]*","",as.character(unlist(c(RawTable[,c(4)])))))
Peptide_Length <- as.character(unlist(c(RawTable[,c(3)])))
Peptide_Mass <- as.character(unlist(c(RawTable[,c(2)])))
Peptide_PTM <- as.character(unlist(c(RawTable[,c(5)])))
Dataset1 <- data.table(Protein_Name,UniProtID,Peptide_Sequence,Peptide_Length,Peptide_PTM,Peptide_Mass)
#DataSet1 <- DS[!duplicated(DS[,c(3)])]

Rawdata <- read.csv(file="~/Documents/Wells Lab/sHLA Results/MCF7/210724_sHLA_MCF7_Senescence_Growing.db.MCF7 Senescence 2.peptides.csv", stringsAsFactors = FALSE)

#Removes everything but the UniProt ID, Peptide Sequence, and Protein Name and counts the Peptide Length
RawTable <- Rawdata[,-c(2,5:15,18:19)]
#RawTable <- Rawdata[-c(1),-c(2,5:16,19:20)]
UniProtID <- gsub("\\|.*","",as.character(unlist(c(RawTable[,c(4)]))))
Peptide_Sequence <- as.character(unlist(c(RawTable[,c(1)])))
#Peptide_Sequence <- str_replace_all(Peptide_Sequence, "[^[:alnum:]]", "")
#Peptide_Sequence <- str_replace_all(Peptide_Sequence, "[[:digit:]]+", "")
Protein_Name <- gsub('^.','',sub("^[^|]*","",as.character(unlist(c(RawTable[,c(4)])))))
Peptide_Length <- as.character(unlist(c(RawTable[,c(3)])))
Peptide_Mass <- as.character(unlist(c(RawTable[,c(2)])))
Peptide_PTM <- as.character(unlist(c(RawTable[,c(5)])))
Dataset2 <- data.table(Protein_Name,UniProtID,Peptide_Sequence,Peptide_Length,Peptide_PTM,Peptide_Mass)
#DataSet1 <- DS[!duplicated(DS[,c(3)])]

Rawdata <- read.csv(file="~/Documents/Wells Lab/sHLA Results/MCF7/210724_sHLA_MCF7_Senescence_Growing.db.MCF7 Senescence 3.peptides.csv", stringsAsFactors = FALSE)

#Removes everything but the UniProt ID, Peptide Sequence, and Protein Name and counts the Peptide Length
RawTable <- Rawdata[,-c(2,5:15,18:19)]
#RawTable <- Rawdata[-c(1),-c(2,5:16,19:20)]
UniProtID <- gsub("\\|.*","",as.character(unlist(c(RawTable[,c(4)]))))
Peptide_Sequence <- as.character(unlist(c(RawTable[,c(1)])))
#Peptide_Sequence <- str_replace_all(Peptide_Sequence, "[^[:alnum:]]", "")
#Peptide_Sequence <- str_replace_all(Peptide_Sequence, "[[:digit:]]+", "")
Protein_Name <- gsub('^.','',sub("^[^|]*","",as.character(unlist(c(RawTable[,c(4)])))))
Peptide_Length <- as.character(unlist(c(RawTable[,c(3)])))
Peptide_Mass <- as.character(unlist(c(RawTable[,c(2)])))
Peptide_PTM <- as.character(unlist(c(RawTable[,c(5)])))
Dataset3 <- data.table(Protein_Name,UniProtID,Peptide_Sequence,Peptide_Length,Peptide_PTM,Peptide_Mass)
#DataSet1 <- DS[!duplicated(DS[,c(3)])]

Rawdata <- read.csv(file="~/Documents/Wells Lab/sHLA Results/MCF7/210724_sHLA_MCF7_Senescence_Growing.db.MCF7 Senescence 4.peptides.csv", stringsAsFactors = FALSE)

#Removes everything but the UniProt ID, Peptide Sequence, and Protein Name and counts the Peptide Length
RawTable <- Rawdata[,-c(2,5:15,18:19)]
#RawTable <- Rawdata[-c(1),-c(2,5:16,19:20)]
UniProtID <- gsub("\\|.*","",as.character(unlist(c(RawTable[,c(4)]))))
Peptide_Sequence <- as.character(unlist(c(RawTable[,c(1)])))
#Peptide_Sequence <- str_replace_all(Peptide_Sequence, "[^[:alnum:]]", "")
#Peptide_Sequence <- str_replace_all(Peptide_Sequence, "[[:digit:]]+", "")
Protein_Name <- gsub('^.','',sub("^[^|]*","",as.character(unlist(c(RawTable[,c(4)])))))
Peptide_Length <- as.character(unlist(c(RawTable[,c(3)])))
Peptide_Mass <- as.character(unlist(c(RawTable[,c(2)])))
Peptide_PTM <- as.character(unlist(c(RawTable[,c(5)])))
Dataset4 <- data.table(Protein_Name,UniProtID,Peptide_Sequence,Peptide_Length,Peptide_PTM,Peptide_Mass)
#DataSet1 <- DS[!duplicated(DS[,c(3)])]
#DS <- data.table(rbind(Dataset1,Dataset2,Dataset3,Dataset4))
DS <- data.table(rbind(Dataset1,Dataset2))
RawTableG <- DS[!duplicated(DS[,c(3)])]


#Will intersect to produce a list of peptides that were found at least 3 times across the 4 senescence datasets
Intersect1a <- merge(Dataset1, Dataset2, by = intersect("Peptide_Sequence", "Peptide_Sequence"),no.dups = FALSE, incomparables = NULL)
Intersect1 <- merge(Intersect1a, Dataset3, by = intersect("Peptide_Sequence", "Peptide_Sequence"),no.dups = FALSE, incomparables = NULL)
Intersect2a <- merge(Dataset1, Dataset3, by = intersect("Peptide_Sequence", "Peptide_Sequence"),no.dups = FALSE, incomparables = NULL)
Intersect2 <- merge(Intersect2a, Dataset4, by = intersect("Peptide_Sequence", "Peptide_Sequence"),no.dups = FALSE, incomparables = NULL)
Intersect3a <- merge(Dataset1, Dataset2, by = intersect("Peptide_Sequence", "Peptide_Sequence"),no.dups = FALSE, incomparables = NULL)
Intersect3 <- merge(Intersect3a, Dataset4, by = intersect("Peptide_Sequence", "Peptide_Sequence"),no.dups = FALSE, incomparables = NULL)
Intersect4a <- merge(Dataset2, Dataset3, by = intersect("Peptide_Sequence", "Peptide_Sequence"),no.dups = FALSE, incomparables = NULL)
Intersect4 <- merge(Intersect4a, Dataset4, by = intersect("Peptide_Sequence", "Peptide_Sequence"),no.dups = FALSE, incomparables = NULL)
Combo <- rbind(Intersect1, Intersect2, Intersect3, Intersect4)
SenescenceCombo <- Combo[!duplicated(Combo[,c(1)])]
write.csv(SenescenceCombo, file="./SenescenceAtLeast3MCF7.csv")

##Will compile a list of all the peptides found in growing
RawdataA <- read.csv(file="~/Documents/Wells Lab/sHLA Results/MCF7/210724_sHLA_MCF7_Senescence_Growing.db.MCF7 Senescence 1.peptides.csv", stringsAsFactors = FALSE)
RawdataB <- read.csv(file="~/Documents/Wells Lab/sHLA Results/MCF7/210724_sHLA_MCF7_Senescence_Growing.db.MCF7 Senescence 2.peptides.csv", stringsAsFactors = FALSE)
RawdataC <- read.csv(file="~/Documents/Wells Lab/sHLA Results/MCF7/210724_sHLA_MCF7_Senescence_Growing.db.MCF7 Senescence 3.peptides.csv", stringsAsFactors = FALSE)
RawdataD <- read.csv(file="~/Documents/Wells Lab/sHLA Results/MCF7/210724_sHLA_MCF7_Senescence_Growing.db.MCF7 Senescence 4.peptides.csv", stringsAsFactors = FALSE)
RawTableA <- RawdataA[,-c(2,5:15,18:19)]
RawTableB <- RawdataB[,-c(2,5:15,18:19)]
RawTableC <- RawdataC[,-c(2,5:15,18:19)]
RawTableD <- RawdataD[,-c(2,5:15,18:19)]
DS <- data.table(rbind(RawTableA,RawTableB,RawTableC,RawTableD))
DS <- data.table(rbind(RawTableA,RawTableB,RawTableC,RawTableD), fill=TRUE)
RawTable <- DS[!duplicated(DS[,c(1)])]
write.csv(RawTable, file="./MCF7GrowingTotal.csv")

#Used Excel to compare peptides in MCF7GrowingTotal.csv and SenescenceAtLeast3MCF7.csv

----------------------------------

#Following script was used for NetMHC data

#Followed same script to parse MS data for B721.221 files

Intersect1 <- merge(merge(Dataset1, Dataset2, by = intersect("Peptide_Sequence", "Peptide_Sequence"),no.dups = FALSE, incomparables = NULL), Dataset3, by = intersect("Peptide_Sequence", "Peptide_Sequence"),no.dups = FALSE, incomparables = NULL)
Intersect2 <- merge(merge(Dataset1, Dataset3, by = intersect("Peptide_Sequence", "Peptide_Sequence"),no.dups = FALSE, incomparables = NULL), Dataset4, by = intersect("Peptide_Sequence", "Peptide_Sequence"),no.dups = FALSE, incomparables = NULL)
Intersect3 <- merge(merge(Dataset1, Dataset2, by = intersect("Peptide_Sequence", "Peptide_Sequence"),no.dups = FALSE, incomparables = NULL), Dataset4, by = intersect("Peptide_Sequence", "Peptide_Sequence"),no.dups = FALSE, incomparables = NULL)
Intersect4 <- merge(merge(Dataset2, Dataset3, by = intersect("Peptide_Sequence", "Peptide_Sequence"),no.dups = FALSE, incomparables = NULL), Dataset4, by = intersect("Peptide_Sequence", "Peptide_Sequence"),no.dups = FALSE, incomparables = NULL)
Combo <- rbind(Intersect1, Intersect2, Intersect3, Intersect4)
Total_Data <- Combo[!duplicated(Combo[,c(1)])]

#Writes FASTA file for 8-mer peptides for NetMHC Analysis
EightMers <-data.table(Total_Data[Peptide_Length.x == "8"])
UNI <- as.list(EightMers$UniProtID.x)
seq <- as.list(EightMers$Peptide_Sequence)
write.fasta(names=UNI, sequences=seq, file.out="Eightmer_prenetmhc2.fasta")

#Writes FASTA file for 9-mer peptides for NetMHC Analysis
EightMers <-data.table(Total_Data[Peptide_Length.x == "9"])
UNI <- as.list(EightMers$UniProtID.x)
seq <- as.list(EightMers$Peptide_Sequence)
write.fasta(names=UNI, sequences=seq, file.out="Ninemer_prenetmhc2.fasta")

#Writes FASTA file for 10-mer peptides for NetMHC Analysis
EightMers <-data.table(Total_Data[Peptide_Length.x == "10"])
UNI <- as.list(EightMers$UniProtID.x)
seq <- as.list(EightMers$Peptide_Sequence)
write.fasta(names=UNI, sequences=seq, file.out="Tenmer_prenetmhc2.fasta")

#Writes FASTA file for 11-mer peptides for NetMHC Analysis
EightMers <-data.table(Total_Data[Peptide_Length.x == "11"])
UNI <- as.list(EightMers$UniProtID.x)
seq <- as.list(EightMers$Peptide_Sequence)
write.fasta(names=UNI, sequences=seq, file.out="Elevenmer_prenetmhc2.fasta")

#Writes FASTA file for 12-mer peptides for NetMHC Analysis
EightMers <-data.table(Total_Data[Peptide_Length.x == "12"])
UNI <- as.list(EightMers$UniProtID.x)
seq <- as.list(EightMers$Peptide_Sequence)
write.fasta(names=UNI, sequences=seq, file.out="Twelvemer_prenetmhc2.fasta")

#You will need to manually run each of these fasta files through the NetMHC4.0 Server and then combine then back into one file

#Removes everything but the UniProt ID, Peptide Sequence, and Protein Name and counts the Peptide Length
NetMHCT8 <- read.csv(file="~/Documents/Wells Lab/sHLA Results/sHLA-B35.1 B721.221 BR1/8MersNetMHC2.csv", stringsAsFactors = FALSE)
NetMHCTable8 <- NetMHCT8[-c(1:2),-c(1)]
NetMHCT9 <- read.csv(file="~/Documents/Wells Lab/sHLA Results/sHLA-B35.1 B721.221 BR1/9MersNetMHC2.csv", stringsAsFactors = FALSE)
NetMHCTable9 <- NetMHCT9[-c(1:2),-c(1)]
NetMHCT10 <- read.csv(file="~/Documents/Wells Lab/sHLA Results/sHLA-B35.1 B721.221 BR1/10MersNetMHC2.csv", stringsAsFactors = FALSE)
NetMHCTable10 <- NetMHCT10[-c(1:2),-c(1)]
NetMHCT11 <- read.csv(file="~/Documents/Wells Lab/sHLA Results/sHLA-B35.1 B721.221 BR1/11MersNetMHC2.csv", stringsAsFactors = FALSE)
NetMHCTable11 <- NetMHCT11[-c(1:2),-c(1)]
NetMHCT12 <- read.csv(file="~/Documents/Wells Lab/sHLA Results/sHLA-B35.1 B721.221 BR1/12MersNetMHC2.csv", stringsAsFactors = FALSE)
NetMHCTable12 <- NetMHCT12[-c(1:2),-c(1)]
NetMHCTable <- rbind(NetMHCTable8,NetMHCTable9,NetMHCTable10,NetMHCTable11,NetMHCTable12)
Peptide <- as.character(unlist(c(NetMHCTable[,c(1)])))
ID <- as.character(unlist(c(NetMHCTable[,c(2)])))
Affinity_nM <- as.numeric(unlist(c(NetMHCTable[,c(3)])))
Rank <- as.numeric(unlist(c(NetMHCTable[,c(4)])))
Number_of_Binders <- as.character(unlist(c(NetMHCTable[,c(7)])))
Prelim_Table <- data.table(Peptide,ID,Affinity_nM,Rank,Number_of_Binders)
#Terminal_Table <- data.table(Prelim_Table[Number_of_Binders == "1"])

#Outputs a file listing all the peptides post-filtering
##Will need to adjust file path for your computer
write.csv(Prelim_Table, file="./NetMHCfilteredpeptides.csv")

#Generates a density plot of the Affinity
countss <- c(unlist(Prelim_Table[,c(3)]))
#dens <- density(countss)
hist(countss, col = "steelblue", frame = FALSE, main = "sHLA-B*35:01 B721.221", xlab = "NetMHC Predicted Affinity (nM)", breaks = 3000, xlim = c(0,2000))
png('printsGreat.png', width = 4, height = 6, units = 'in', res = 300)
#polygon(dens, col = "steelblue")

#Generates a density plot of the Affinity
countss <- c(unlist(Prelim_Table[,c(3)])) 
dens <- density(countss)
plot(dens, col = "steelblue", frame = FALSE, main = "sHLA-B*35:01 B721.221", xlab = "NetMHC Predicted Affinity (nM)", xlim = c(0,2000))
polygon(dens, col = "steelblue")

#Generates a density plot of the %Rank
countss <- c(unlist(Prelim_Table[,c(4)])) 
dens <- density(countss)
plot(dens, col = "steelblue", frame = FALSE, main = "sHLA-B*35:01 B721.221", xlab = "NetMHC %Rank", xlim = c(0,3), ylim = c(0,1))
polygon(dens, col = "steelblue")