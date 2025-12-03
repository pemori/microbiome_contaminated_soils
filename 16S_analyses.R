#Scripts for the article "Microbial refuge under shrubs limits soil Cu contamination effects"
#Submitted to Soil Biology and Biochemistry
#
#December 2nd, 2025

# Inquiries to Dr. Pedro Mondaca > pedro.mondaca@pucv.cl


setwd("C:\\Users\\pedro\\Research\\microbiome_contaminated_soils\\Github")

#Load libraries
{
library(ggplot2)  
library(readxl)  
library(vegan)
library(dplyr)
library(phyloseq)
library(microeco)
library(microbiome)
library(magrittr)
library(ggplot2)
library(data.table)
library(ggstatsplot)
library(ggsignif)
library(rcompanion)
library(userfriendlyscience)
library(patchwork)
library(pairwiseAdonis)
library(pheatmap)
library(ALDEx2)
library(igraph)
library(rgexf)
library(WGCNA)
library(NetCoMi)
library(readxl)  
library(vegan)
library(dplyr)
library(tidyr)
library(phyloseq)
library(corrplot)
library(microeco)
library(microbiome)
library(magrittr)
library(ggplot2)
library(data.table)
library(ggstatsplot)
library(ggsignif)
library(rcompanion)
library(patchwork)
library(pairwiseAdonis)
library(fantaxtic)
library(pheatmap)
library(ALDEx2)
library(igraph)
library(rgexf)
library(WGCNA)
library(NetCoMi)
library(reshape2)
library(MicroNiche)
library(userfriendlyscience)
}


#Load data ####
{
abund_seq <- read.csv("ASV_abundance_Table_16S.csv",sep=',',header = FALSE)
seqs<-abund_seq[1,]
abund <- read.csv("ASV_abundance_Table_16S.csv",sep=',',header = TRUE)
rownames(abund)<-abund[,1]
abund<-abund[-1]
dim(abund)
abund[1:5, 15820:15823]
#otu_table ready

tax<-read.csv("ASV_tax_assignments_16s.txt",header=FALSE)
head(tax)
library(tidyr)
tax_sep <-tax %>% separate(V1,c("seqs", "tax"), sep = "k_") %>%
  dplyr::mutate(seqs = stringr::str_remove(seqs, pattern = "\t"))
tax_sep$tax <-paste("k",tax_sep$tax, sep="")
tax_sep <-tax_sep %>% separate(tax,c("Kingdom","Phylum","Class","Order","Family","Genus","Specie"), sep = ";")
head(tax_sep)
row.names(tax_sep)<-tax_sep[,1]
tax_num<-tax_sep[-1]
dim(tax_num)
tax_num[15345:15349, 1:7]
#tax_table almost ready

#match the rows
length(colnames(abund))
genomic_idx <- match(colnames(abund), rownames(tax_num))
length(genomic_idx)
genomic_idx
tax_ord <-tax_num[genomic_idx,]
tax_ord
rownames(tax_ord) <- colnames(abund)
dim(tax_ord)
#tax_table ready

library(dplyr)
sample<-factor(c("Bulk Soil","Bulk Soil","Bulk Soil","Soil surrounding root","Soil surrounding root","Soil surrounding root"))
contamination<-factor(c("High Cont","Mid Cont","Uncont","High Cont","Mid Cont","Uncont"))
sampxcont<-factor(c("BS HC","BS MC","BS UC","SSR HC","SSR MC","SSR UC"))
mapping<-data.frame(sample,contamination,sampxcont)
mapping <-mapping %>% slice(rep(1:n(),each=10))
row.names(mapping) <- row.names(abund)
mapping
#sample_data ready

library(phyloseq)
OTU = otu_table(as.matrix(abund), taxa_are_rows = FALSE) # create the occurrence table object in phyloseq format
SAM = sample_data(mapping, errorIfNULL = TRUE) # create the sample metadata object in phyloseq format
TAX = tax_table(as.matrix(tax_ord)) # create the observation metadata object (OTU taxonomy) in phyloseq format
data_phylo0 <- phyloseq(OTU, TAX, SAM) # create the phyloseq object including occurrence table data and sample/observation metadata
data_phylo0 # print information about the phyloseq object

filterPhyla="k_Bacteria"
data_phylo = subset_taxa(data_phylo0, Kingdom %in% filterPhyla) 
#phyloseq file ready

#data processing
nb_samples <- dim(abund)[1] # nb of rows, here samples
nb_samples
nb_var <- dim(abund)[2] # number of columns, here variables (OTU)
nb_var
sum(abund == 0)
sum(abund == 0) / (nb_var * nb_samples) * 100 #96% of the data are zeros

hist(as.matrix(abund),
     right = FALSE, 
     las = 1, 
     xlab = "Occurrence value", 
     ylab = "Frequency", 
     main = "Occurrence frequency")
min(colSums(abund))

#In order to check how the different OTU/ASV are shared between samples, plot the number of non-zero values for each OTU
non_zero <- 0*1:nb_var
for (i in 1:nb_var){
  non_zero[i]<-sum(abund[,i] != 0)
}
plot(non_zero, xlab = "ASV's", ylab = "Frequency", main = "Number of non zero values", las = 1)
rarecurve(abund, step = 100, cex = 0.75, las = 1) 

#Library size
sum_seq <- rowSums(abund)
plot(sum_seq, ylim=c(0,80000), main=c("Number of counts per sample"), xlab=c("Samples"))
sum_seq
min(sum_seq)
max(sum_seq)

# filter the OTU data using filter_taxa function included in phyloseq package
data_phylo_filt = filter_taxa(data_phylo, function(x) sum(x > 2) > (0.033 * length(x)), TRUE) # 2/n= 2/60 = 0.033
data_otu_filt = data.frame(otu_table(data_phylo_filt)) 
dim(data_otu_filt)
dim(abund)
1-(4916/15823) # so I removed 68.9 % of the ASVs
sum(data_otu_filt)
sum(abund)
((sum(abund)-sum(data_otu_filt))/sum(abund))*100 # But only 21% of the total reads. 
#That is, I removed those extremely rare ASVs and present in just few samples

#Rarefaction
otu_mat <- as(otu_table(data_phylo_filt), "matrix")
if(taxa_are_rows(data_phylo_filt)) {
  otu_mat <- t(otu_mat)
}
rarecurve(otu_mat, step = 100, sample = min(rowSums(otu_mat)), 
          col = "steelblue", cex = 0.6, label = FALSE,
          ylab = "Observed ASVs", xlab = "Sequencing Depth")
sample_sums(data_phylo_filt)
summary(sample_sums(data_phylo_filt))  # for min, max, median, mean, etc.

#rarefacción
set.seed(2) # set seed for analysis reproducibility
data_phylo_filt_rar <- phyloseq::rarefy_even_depth(
  data_phylo_filt,
  sample.size = 11892,    
  rngseed = 2,          
  replace = FALSE,
  trimOTUs = TRUE,
  verbose = TRUE
)

data_phylo_filt_rar <- subset_taxa(data_phylo_filt_rar, Genus != "g__soil")
data_phylo_filt_rar <- subset_taxa(data_phylo_filt_rar, Family != "f__Mitochondria")
data_phylo_filt_rar <- subset_taxa(data_phylo_filt_rar, Class != "c__Chloroplast")
data_otu_filt_rar = data.frame(otu_table(data_phylo_filt_rar)) 
sum_seq_rar <- rowSums(data_otu_filt_rar)
sum_seq_rar #11892 reads

data_phylo_filt_rar_gen <- tax_glom(data_phylo_filt_rar, taxrank = "Genus")
genus_values <- tax_table(data_phylo_filt_rar_gen)[, "Genus"]
family_values <- tax_table(data_phylo_filt_rar_gen)[, "Family"]
order_values <- tax_table(data_phylo_filt_rar_gen)[, "Order"]
class_values <- tax_table(data_phylo_filt_rar_gen)[, "Class"]
phylum_values <- tax_table(data_phylo_filt_rar_gen)[, "Phylum"]
bacteria_values <- tax_table(data_phylo_filt_rar_gen)[, "Kingdom"]

# Replace 'g__NA' in 'genus' with corresponding 'family' values
genus_values[genus_values == "g__NA"] <- family_values[genus_values == "g__NA"]
genus_values[genus_values == "f__NA"] <- order_values[genus_values == "f__NA"]
genus_values[genus_values == "o__NA"] <- class_values[genus_values == "o__NA"]
genus_values[genus_values == "c__NA"] <- phylum_values[genus_values == "c__NA"]
genus_values[genus_values == "p__NA"] <- bacteria_values[genus_values == "p__NA"]
# Assign the updated 'genus' column back to the tax_table in the phyloseq object
tax_table(data_phylo_filt_rar_gen)[, "Genus"] <- genus_values
tax_table <- data.frame(tax_table(data_phylo_filt_rar_gen))

# Check the tax_table to verify changes
tax_table
tax_table_filt_rar_gen <- tax_table[,-7]
tax_table_filt_rar_gen
genus <- as.character(tax_table_filt_rar_gen$Genus)
genus
bact_positions <- which(genus == "g__bacterium")
conex_positions<- which(genus == "g__Conexibacter")
conex_positions
genus[bact_positions] <- paste0("g__bacterium_", seq_along(bact_positions))
genus[conex_positions] <- paste0("g__Conexibacter_", seq_along(conex_positions))
genus <- sub("^g__", "", genus)
genus <- sub("^f__", "", genus)
genus <- sub("^o__", "", genus)
genus <- sub("^c__", "", genus)
genus <- sub("^p__", "", genus)
genus
all_duplicates <- duplicated(tax_table$Genus, fromLast = TRUE)
print(all_duplicates)

Acidobac_positions <- which(genus == "Acidobacteria")
Actinobac_positions<- which(genus == "Actinobacteria")
genus[Acidobac_positions] <- paste0("Acidobacteria_", seq_along(bact_positions))
genus[Actinobac_positions] <- paste0("Actinobacteria_", seq_along(conex_positions))
genus

row.names(tax_table_filt_rar_gen)<-genus
tax_table_filt_rar_gen

data_otu_filt_rar_gen = data.frame(otu_table(data_phylo_filt_rar_gen)) 
data_otu_filt_rar_gen
colnames(data_otu_filt_rar_gen)<-genus
data_otu_filt_rar_gen

OTU_g = otu_table(as.matrix(data_otu_filt_rar_gen), taxa_are_rows = FALSE) # create the occurrence table object in phyloseq format
SAM_g = sample_data(mapping, errorIfNULL = TRUE) # create the sample metadata object in phyloseq format
TAX_g = tax_table(as.matrix(tax_table_filt_rar_gen)) # create the observation metadata object (OTU taxonomy) in phyloseq format
data_phylo_g <- phyloseq(OTU_g, TAX_g, SAM_g) # create the phyloseq object including occurrence table data and sample/observation metadata
data_phylo_g # print information about the phyloseq object

data_phylo_b <- data_phylo_filt_rar
data_phylo_bg <- data_phylo_g

phylo_BS_UC <- subset_samples(data_phylo_filt_rar, sampxcont == "BS UC")
phylo_BS_MC <- subset_samples(data_phylo_filt_rar, sampxcont == "BS MC")
phylo_BS_HC <- subset_samples(data_phylo_filt_rar, sampxcont == "BS HC")
phylo_SSR_UC <- subset_samples(data_phylo_filt_rar, sampxcont == "SSR UC")
phylo_SSR_MC <- subset_samples(data_phylo_filt_rar, sampxcont == "SSR MC")
phylo_SSR_HC <- subset_samples(data_phylo_filt_rar, sampxcont == "SSR HC")
phylo_BS <- subset_samples(data_phylo_filt_rar, sample == "Bulk Soil")
phylo_SSR <- subset_samples(data_phylo_filt_rar, sample == "Soil surrounding root")
otu_BS_UC = data.frame(otu_table(phylo_BS_UC))
otu_BS_MC = data.frame(otu_table(phylo_BS_MC))
otu_BS_HC = data.frame(otu_table(phylo_BS_HC))
otu_SSR_UC = data.frame(otu_table(phylo_SSR_UC))
otu_SSR_MC = data.frame(otu_table(phylo_SSR_MC))
otu_SSR_HC = data.frame(otu_table(phylo_SSR_HC))
otu_BS = data.frame(otu_table(phylo_BS))
otu_SSR = data.frame(otu_table(phylo_SSR))      

phylo_BS_UC_g <- subset_samples(data_phylo_g, sampxcont == "BS UC")
phylo_BS_MC_g <- subset_samples(data_phylo_g, sampxcont == "BS MC")
phylo_BS_HC_g <- subset_samples(data_phylo_g, sampxcont == "BS HC")
phylo_SSR_UC_g <- subset_samples(data_phylo_g, sampxcont == "SSR UC")
phylo_SSR_MC_g <- subset_samples(data_phylo_g, sampxcont == "SSR MC")
phylo_SSR_HC_g <- subset_samples(data_phylo_g, sampxcont == "SSR HC")
phylo_BS_g <- subset_samples(data_phylo_g, sample == "Bulk Soil")
phylo_SSR_g <- subset_samples(data_phylo_g, sample == "Soil surrounding root")
otu_BS_UC_g = data.frame(otu_table(phylo_BS_UC_g))
otu_BS_MC_g = data.frame(otu_table(phylo_BS_MC_g))
otu_BS_HC_g = data.frame(otu_table(phylo_BS_HC_g))
otu_SSR_UC_g = data.frame(otu_table(phylo_SSR_UC_g))
otu_SSR_MC_g = data.frame(otu_table(phylo_SSR_MC_g))
otu_SSR_HC_g = data.frame(otu_table(phylo_SSR_HC_g))
otu_BS_g = data.frame(otu_table(phylo_BS_g))
otu_SSR_g = data.frame(otu_table(phylo_SSR_g))

#Now I will prepare the data for subsequent niche breadth analyses
tax_table(data_phylo_filt_rar)
data_phylo_filt_rar_fam_b <- tax_glom(data_phylo_filt_rar, taxrank = "Family") #agglomerated but not unique names
# # Access the taxonomic table and get the Family column
Family_values <- tax_table(data_phylo_filt_rar_fam_b)[, "Family"]
order_values <- tax_table(data_phylo_filt_rar_fam_b)[, "Order"]
class_values <- tax_table(data_phylo_filt_rar_fam_b)[, "Class"]
phylum_values <- tax_table(data_phylo_filt_rar_fam_b)[, "Phylum"]
bacteria_values <- tax_table(data_phylo_filt_rar_fam_b)[, "Kingdom"]
Family_values
# Replace 'g__NA' in 'Family' with corresponding 'family' values
Family_values[Family_values == "f__NA"] <- order_values[Family_values == "f__NA"]
Family_values[Family_values == "o__NA"] <- class_values[Family_values == "o__NA"]
Family_values[Family_values == "c__NA"] <- phylum_values[Family_values == "c__NA"]
Family_values[Family_values == "p__NA"] <- bacteria_values[Family_values == "p__NA"]
Family_values

# Assign the updated 'Family' column back to the tax_table in the phyloseq object
tax_table(data_phylo_filt_rar_fam_b)[, "Family"] <- Family_values
tax_table <- data.frame(tax_table(data_phylo_filt_rar_fam_b))

# Check the tax_table to verify changes
tax_table
tax_table_filt_rar_fam_b <- tax_table[,-c(6:7)]
tax_table_filt_rar_fam_b

Family <- as.character(tax_table_filt_rar_fam_b$Family)
Family

data_phylo_filt_rar_fam_b
tax_fam_b <- as.data.frame(tax_table(data_phylo_filt_rar_fam_b))
dim(tax_fam_b)
row.names(tax_table_filt_rar_fam_b)<-Family
tax_table_filt_rar_fam_b

data_otu_filt_rar_fam_b = data.frame(otu_table(data_phylo_filt_rar_fam_b)) 
data_otu_filt_rar_fam_b
tax_fam_b
colnames(data_otu_filt_rar_fam_b)<-tax_fam_b$Family
data_otu_filt_rar_fam_b

# tax_table(data_phylo_filt_rar_fam_b)

#new phyloseq
OTU = otu_table(as.matrix(data_otu_filt_rar_fam_b), taxa_are_rows = FALSE) # create the occurrence table object in phyloseq format
SAM = sample_data(mapping, errorIfNULL = TRUE) # create the sample metadata object in phyloseq format
TAX = tax_table(as.matrix(tax_table_filt_rar_fam_b)) # create the observation metadata object (OTU taxonomy) in phyloseq format
phylo_bac_f <- phyloseq(OTU, TAX, SAM)
phylo_bac_f


# Filter by category
Uncont_fam_b <- data_otu_filt_rar_fam_b[grepl("Uncont", rownames(data_otu_filt_rar_fam_b)), ]
SC_fam_b <- data_otu_filt_rar_fam_b[grepl("SC", rownames(data_otu_filt_rar_fam_b)), ]
Cont_fam_b <- data_otu_filt_rar_fam_b[grepl("Cont", rownames(data_otu_filt_rar_fam_b)), ]
SS_fam_b <- data_otu_filt_rar_fam_b[!grepl("SSR", rownames(data_otu_filt_rar_fam_b)), ]
SSR_fam_b <- data_otu_filt_rar_fam_b[grepl("SSR", rownames(data_otu_filt_rar_fam_b)), ]
SSR_fam_b

}

#Alpha diversity####
{
#View(data_otu_filt_rar)  
data_richness <- estimateR(data_otu_filt_rar)                                            # calculate richness and Chao1 using vegan package
R_family <- specnumber(data_otu_filt_rar)
data_evenness <- vegan::diversity(data_otu_filt_rar) / log(specnumber(data_otu_filt_rar))                # calculate evenness index using vegan package
data_shannon <- vegan::diversity(data_otu_filt_rar, index = "shannon")                          # calculate Shannon index using vegan package
data_richness2<-t(data_richness)
data_richness2 
tab_diver<-cbind(R_family,data_evenness,data_shannon,data_richness2) 
tab_diver 
dim(mapping)
dim(tab_diver)
tab_diver2<-cbind(mapping, tab_diver)
tab_diver2<-as.data.frame(tab_diver2)
str(tab_diver2)
#library(xlsx)
#write.xlsx(tab_diver2,"tab_diver2.xlsx",sheetName = "diver")

BGQ <- read_excel("BGQ3.xlsx", 
                   sheet = "All", col_types = c("text", "text",
                                                "text", "text", "text", "numeric","numeric", 
                                                "numeric", "numeric", "numeric", 
                                                "numeric", "numeric", "numeric", 
                                                "numeric", "numeric"))
dim(tab_diver2)
dim(BGQ)

tab_diver2_avg <- tab_diver2 %>%
  group_by(sampxcont) %>%
  summarise(across(c(R_family, data_evenness, data_shannon, 
                     S.obs, S.chao1, se.chao1, S.ACE, se.ACE),
                   mean, na.rm = TRUE))

tab_diver2_avg

# Reorder the factor levels
order_contamination <- c("Uncont", "Mid Cont", "High Cont")
colors_contamination <- c("darkgreen","gold", "red")
new_labels <- c("BS", "SSR")

#diversity plots
s_box<-ggplot(tab_diver2, aes(x = sample, y = S.obs, fill = factor(contamination, levels = order_contamination))) +
  geom_boxplot() +
  labs(x = " ", y = "Richness") +
  scale_fill_manual(values = colors_contamination) +
  scale_y_continuous(limits = c(0, 700))+
  theme_bw() + guides(fill = "none") +
  theme(axis.text.x = element_text(size = 12, color = "black", family = "Arial"),
        axis.text.y = element_text(size = 11, color = "black", family = "Arial"),
        axis.title.y = element_text(size = 14),
        panel.grid.major = element_line(color = "white", linetype = "dashed"))+
  scale_x_discrete(labels=new_labels)
s_box

chao_box<-ggplot(tab_diver2, aes(x = sample, y = S.chao1, fill = factor(contamination, levels = order_contamination))) +
  geom_boxplot() +
  labs(x = " ", y = "Chao1 Index") +
  scale_fill_manual(values = colors_contamination) +
  #scale_y_continuous(limits = c(3, 7))+
  theme_bw() + guides(fill = "none") +
  theme(axis.text.x = element_text(size = 12, color = "black", family = "Arial"),
        axis.text.y = element_text(size = 11, color = "black", family = "Arial"),
        axis.title.y = element_text(size = 14),
        panel.grid.major = element_line(color = "white", linetype = "dashed"))+
  scale_x_discrete(labels=new_labels)
chao_box

shan_box<-ggplot(tab_diver2, aes(x = sample, y = data_shannon, fill = factor(contamination, levels = order_contamination))) +
  geom_boxplot() +
  labs(x = " ", y = "Shannon entropy") +
  scale_fill_manual(values = colors_contamination) +
  scale_y_continuous(limits = c(3, 7))+
  theme_bw() + guides(fill = "none") +
  theme(axis.text.x = element_text(size = 12, color = "black", family = "Arial"),
        axis.text.y = element_text(size = 11, color = "black", family = "Arial"),
        axis.title.y = element_text(size = 14),
        panel.grid.major = element_line(color = "white", linetype = "dashed"))+
  scale_x_discrete(labels=new_labels)
shan_box

even_box<-ggplot(tab_diver2, aes(x = sample, y = data_evenness, fill = factor(contamination, levels = order_contamination))) +
  geom_boxplot() +
  labs(x = " ", y = "Evenness") +
  scale_fill_manual(values = colors_contamination) +
 # scale_y_continuous(limits = c(0.75, 1))+
  theme_bw() + guides(fill = "none") +
  theme(axis.text.x = element_text(size = 12, color = "black", family = "Arial"),
        axis.text.y = element_text(size = 11, color = "black", family = "Arial"),
        axis.title.y = element_text(size = 14),
        panel.grid.major = element_line(color = "white", linetype = "dashed"))+
  scale_x_discrete(labels=new_labels)
even_box

library(patchwork)
(s_box | shan_box | even_box) #printed at 950x360

#let's separe df for two-way ANOVA
#first, comparing contamination levels
str(tab_diver2)
BS<-tab_diver2 %>% filter(sample=="Bulk Soil")
SSR<-tab_diver2 %>% filter(sample=="Soil surrounding root")

#BS
# ANOVA-like test (Welch ANOVA)
oneway.test(S.obs ~ contamination, data = BS, var.equal = FALSE)
# Games-Howell post hoc test
games_howell_result <- posthocTGH(BS$S.obs, BS$contamination, method="games-howell")
print(games_howell_result)

# ANOVA-like test (Welch ANOVA)
oneway.test(S.chao1 ~ contamination, data = BS, var.equal = FALSE)
# Games-Howell post hoc test
games_howell_result <- posthocTGH(BS$S.obs, BS$contamination, method="games-howell")
print(games_howell_result)

# ANOVA-like test (Welch ANOVA)
oneway.test(data_shannon ~ contamination, data = BS, var.equal = FALSE)
# Games-Howell post hoc test
games_howell_result <- posthocTGH(BS$data_shannon, BS$contamination, method="games-howell")
print(games_howell_result)

# ANOVA-like test (Welch ANOVA)
oneway.test(data_evenness ~ contamination, data = BS, var.equal = FALSE)
# Games-Howell post hoc test
games_howell_result <- posthocTGH(BS$data_evenness, BS$contamination, method="games-howell")
print(games_howell_result)

#SSR
# ANOVA-like test (Welch ANOVA)
oneway.test(S.obs ~ contamination, data = SSR, var.equal = FALSE)
# Games-Howell post hoc test
games_howell_result <- posthocTGH(SSR$S.obs, SSR$contamination, method="games-howell")
print(games_howell_result)

# ANOVA-like test (Welch ANOVA)
oneway.test(data_shannon ~ contamination, data = SSR, var.equal = FALSE)
# Games-Howell post hoc test
games_howell_result <- posthocTGH(SSR$data_shannon, SSR$contamination, method="games-howell")
print(games_howell_result)

# ANOVA-like test (Welch ANOVA)
oneway.test(data_evenness ~ contamination, data = SSR, var.equal = FALSE)
# Games-Howell post hoc test
games_howell_result <- posthocTGH(SSR$data_evenness, SSR$contamination, method="games-howell")
print(games_howell_result)

#now, by comparing soil compartments: BS y SSR
str(tab_diver2)
UC<-tab_diver2 %>% filter(contamination=="Uncont")
MC<-tab_diver2 %>% filter(contamination=="Mid Cont")
HC<-tab_diver2 %>% filter(contamination=="High Cont")

# ANOVA-like test (Welch ANOVA)
oneway.test(S.obs ~ sample, data = UC, var.equal = FALSE)
# Games-Howell post hoc test
games_howell_result <- posthocTGH(UC$S.obs, UC$sample, method="games-howell")
print(games_howell_result)

# ANOVA-like test (Welch ANOVA)
oneway.test(data_shannon ~ sample, data = UC, var.equal = FALSE)
# Games-Howell post hoc test
games_howell_result <- posthocTGH(UC$data_shannon, UC$sample, method="games-howell")
print(games_howell_result)

# ANOVA-like test (Welch ANOVA)
oneway.test(data_evenness ~ sample, data = UC, var.equal = FALSE)
# Games-Howell post hoc test
games_howell_result <- posthocTGH(UC$data_evenness, UC$sample, method="games-howell")
print(games_howell_result)

# ANOVA-like test (Welch ANOVA)
oneway.test(S.obs ~ sample, data = MC, var.equal = FALSE)
# Games-Howell post hoc test
games_howell_result <- posthocTGH(MC$S.obs, MC$sample, method="games-howell")
print(games_howell_result)

# ANOVA-like test (Welch ANOVA)
oneway.test(data_shannon ~ sample, data = MC, var.equal = FALSE)
# Games-Howell post hoc test
games_howell_result <- posthocTGH(MC$data_shannon, MC$sample, method="games-howell")
print(games_howell_result)

# ANOVA-like test (Welch ANOVA)
oneway.test(data_evenness ~ sample, data = MC, var.equal = FALSE)
# Games-Howell post hoc test
games_howell_result <- posthocTGH(MC$data_evenness, MC$sample, method="games-howell")
print(games_howell_result)

# ANOVA-like test (Welch ANOVA)
oneway.test(S.obs ~ sample, data = HC, var.equal = FALSE)
# Games-Howell post hoc test
games_howell_result <- posthocTGH(HC$S.obs, HC$sample, method="games-howell")
print(games_howell_result)

# ANOVA-like test (Welch ANOVA)
oneway.test(data_shannon ~ sample, data = HC, var.equal = FALSE)
# Games-Howell post hoc test
games_howell_result <- posthocTGH(HC$data_shannon, HC$sample, method="games-howell")
print(games_howell_result)

# ANOVA-like test (Welch ANOVA)
oneway.test(data_evenness ~ sample, data = HC, var.equal = FALSE)
# Games-Howell post hoc test
games_howell_result <- posthocTGH(HC$data_evenness, HC$sample, method="games-howell")
print(games_howell_result)

}

#beta-diversity####
{
dist_jaccard <- as.matrix(vegdist(data_otu_filt_rar, method = "jaccard")) 
dist_bc <- as.matrix(vegdist(data_otu_filt_rar, method = "bray")) 

# calculate PCOA using Phyloseq package
pcoa_jaccard = ordinate(data_phylo_filt_rar, "PCoA", "jaccard")
pcoa_bc = ordinate(data_phylo_filt_rar, "PCoA", "bray") 

#plot
plot_ordination(data_phylo_filt_rar, pcoa_bc, color = "contamination", shape = "sample") + 
  geom_point(size = 3) +
  scale_color_manual(values = colors_contamination, breaks = order_contamination) +
  #scale_shape_manual(values = c(1, 2, 3)) +  # Change to suitable shapes
  theme_bw() +
  theme(panel.grid.major = element_line(color = "white", linetype = "dashed"))
#printed at 614 * 422

# PERMANOVA (Permutational Multivariate Analysis of Variance)
# Permanova test using the vegan package 
adonis2(data_otu_filt_rar~contamination, data=mapping, permutations=9999, method="bray") #sign
adonis2(data_otu_filt_rar~sample, data=mapping, permutations=9999, method="bray") #sign

dist_matrix <- vegdist(data_otu_filt_rar, method = "bray")
#Permutational Homogeneity of Multivariate Dispersions
permdisp_result <- betadisper(dist_matrix, mapping$contamination)
permutest(permdisp_result, pairwise = TRUE, permutations = 999)
permdisp_result <- betadisper(dist_matrix, mapping$sample)
permutest(permdisp_result, pairwise = TRUE, permutations = 999)

#homogeneity of multivariate dispersion are not equal among the level of the factors
#Thus, post-hoc difference can also be attributed to dispersion 
#https://besjournals.onlinelibrary.wiley.com/doi/10.1111/j.2041-210X.2011.00127.x

# for posthoc analysis, I will separate the df
BS<-data_otu_filt_rar[1:30,]
BS_map<-mapping[1:30,]
SSR<-data_otu_filt_rar[31:60,]
SSR_map<-mapping[31:60,]
BS_map
HC<-data_otu_filt_rar[c(1:10,31:40),]
HC_map<-mapping[c(1:10,31:40),]
MC<-data_otu_filt_rar[c(11:20,41:50),]
MC_map<-mapping[c(11:20,41:50),]
UC<-data_otu_filt_rar[c(21:30,51:60),]
UC_map<-mapping[c(21:30,51:60),]

library(pairwiseAdonis)
pairwise.adonis2(BS ~ contamination, data = BS_map)
pairwise.adonis2(SSR ~ contamination, data = SSR_map)
pairwise.adonis2(UC ~ sample, data = UC_map)
pairwise.adonis2(MC ~ sample, data = MC_map)
pairwise.adonis2(HC ~ sample, data = HC_map)

}

# Compositional and functional analyses ####
# data_phylo_filt_rar = phyloseq file filtrado y rarefied
# Composition analysis
{
data_phylo_filt_rar_fam_b
tax_data = data.frame(tax_table(data_phylo_filt_rar_fam_b)) # create a separated file
otu_data = (data.frame(otu_table(data_phylo_filt_rar_fam_b)))
sam_data = data.frame(sample_data(data_phylo_filt_rar_fam_b))
otu_data_t <- as.data.frame(t(otu_data))
sam_data<-as.data.frame(sam_data)
rownames(sam_data)<-colnames(otu_data_t)
sam_data
dim(otu_data_t)
dataset <- microtable$new(sample_table = sam_data, otu_table = otu_data_t, tax_table = tax_data)
dataset
dataset$tax_table %<>% base::subset(Kingdom == "k_Bacteria" )
dataset$tidy_dataset()
#all is bacteria

dataset$cal_abund() 
dataset$taxa_abund$Phylum[1:5, 1:5]
dataset$taxa_abund$Family[1:5, 1:5]
dataset$cal_alphadiv(PD = FALSE) #PD=TRUE para phylogenetic analyses
dataset$alpha_diversity
dataset$cal_betadiv(unifrac = FALSE)
dataset$beta_diversity
dataset$tax_table
dataset$sample_table$sampxcont <- factor(dataset$sample_table$sampxcont, levels = c("BS UC", "BS MC", "BS HC", "SSR UC", "SSR MC", "SSR HC"))
dataset$sample_table$sample <- factor(dataset$sample_table$sample, levels = c("Bulk Soil", "Soil surrounding root"))
dataset$sample_table$contamination <- factor(dataset$sample_table$contamination, levels = c("Uncont", "Mid Cont","High Cont"))
dataset

BS_mt <- clone(dataset)
BS_mt$sample_table <- subset(BS_mt$sample_table, sample == "Bulk Soil")
BS_mt$otu_table <- BS_mt$otu_table[,1:30]
BS_mt$tax_table
BS_mt$otu_table
BS_mt$tidy_dataset()
BS_mt

SSR_mt <- clone(dataset)
SSR_mt$sample_table <- subset(SSR_mt$sample_table, sample == "Soil surrounding root")
SSR_mt$otu_table <- SSR_mt$otu_table[,31:60]
SSR_mt$tidy_dataset()
SSR_mt

UC_mt <- clone(dataset)
UC_mt$sample_table <- subset(UC_mt$sample_table, contamination == "Uncont")
UC_mt$otu_table <- UC_mt$otu_table[,c(21:30,51:60)]
UC_mt$tidy_dataset()
UC_mt

MC_mt <- clone(dataset)
MC_mt$sample_table <- subset(MC_mt$sample_table, contamination == "Mid Cont")
MC_mt$otu_table <- MC_mt$otu_table[,c(11:20,41:50)]
MC_mt$tidy_dataset()
MC_mt

HC_mt <- clone(dataset)
HC_mt$sample_table <- subset(HC_mt$sample_table, contamination == "High Cont")
HC_mt$otu_table <- HC_mt$otu_table[,c(1:10,31:40)]
HC_mt$tidy_dataset()
HC_mt

#composition-based class

# show 30 taxa at Family level
t1 <- trans_abund$new(dataset = dataset, taxrank = "Family", ntaxa = 30)
t1$plot_heatmap(facet = c("sample","contamination"), xtext_keep = FALSE, withmargin = FALSE) #836*680


t1 <- trans_abund$new(dataset = dataset, taxrank = "Phylum", ntaxa = 8)
library(ggh4x)
t1$plot_bar(others_color = "grey70", facet = c("sample","contamination"), xtext_keep = FALSE, legend_text_italic = FALSE)
# printed at 1600x800 rel_abund

# The groupmean parameter can be used to obtain the group-mean barplot.
t1 <- trans_abund$new(dataset = dataset, taxrank = "Phylum", ntaxa = 8, groupmean = "sampxcont")
g1 <- t1$plot_bar(others_color = "grey70", legend_text_italic = FALSE)
g1 + theme_classic() + theme(axis.title.y = element_text(size = 18))
#printed at 578*426 mean_rel_abund

#ALDEx2
data_BS <- clone(dataset)
data_BS$sample_table <- subset(data_BS$sample_table, sample == "Bulk Soil")
data_BS$tidy_dataset()
data_BS$sample_table$sample

data_SSR <- clone(dataset)
data_SSR$sample_table <- subset(data_SSR$sample_table, sample == "Soil surrounding root")
data_SSR$tidy_dataset()
data_SSR

data_UC <- clone(dataset)
data_UC$sample_table <- subset(data_UC$sample_table, contamination == "Uncont")
data_UC$tidy_dataset()
data_UC$sample_table$sample
data_UC

data_MC <- clone(dataset)
data_MC$sample_table <- subset(data_MC$sample_table, contamination == "Mid Cont")
data_MC$tidy_dataset()
data_MC$sample_table$sample
data_MC

data_HC <- clone(dataset)
data_HC$sample_table <- subset(data_HC$sample_table, contamination == "High Cont")
data_HC$tidy_dataset()
data_HC$sample_table$sample
data_HC

# ALDEx2_kw for UC
###
library(reshape2)
library(dplyr)
library(ggplot2)

t1 <- trans_diff$new(dataset = data_UC, method = "ALDEx2_kw", group = "sample", taxa_level = "Family")
t1$plot_diff_abund(use_number = 1:11, group_order = c("Bulk Soil","Soil surrounding root"))

data_UC$cal_abund()
abund_family <- data_UC$taxa_abund$Family
abund_family
abund_long <- melt(as.matrix(abund_family))
colnames(abund_long) <- c("Family", "SampleID", "Abundance")
abund_long$Family <- sub(".*f__", "", abund_long$Family)

df_diff <- data.frame(
  Taxa = t1$res_diff$Taxa,
  P.adj = t1$res_diff$P.adj,
  stringsAsFactors = FALSE
)
df_diff <- df_diff[!grepl("f__NA", df_diff$Taxa), ]
df_diff$Family <- sub(".*f__", "", df_diff$Taxa)

top_families <- df_diff %>%
  arrange(P.adj) %>%
  slice(1:10) %>%
  pull(Family)


df_diff %>%
  filter(Family %in% top_families) %>%
  arrange(P.adj)
abund_long$Family <- factor(abund_long$Family, levels = rev(top_families))

abund_plot <- abund_long %>%
  filter(Family %in% top_families)

meta <- data_UC$sample_table
meta$SampleID <- rownames(meta)
meta
abund_plot <- abund_plot %>%
  left_join(meta, by = "SampleID")

abund_summary <- abund_plot %>%
  group_by(Family, sample) %>%
  summarise(
    mean_abund = mean(Abundance, na.rm = TRUE),
    se_abund = sd(Abundance, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

abund_summary$Family <- factor(abund_summary$Family, levels = levels(abund_plot$Family))
abund_summary$sample <- factor(abund_summary$sample, levels = c("Soil surrounding root","Bulk Soil"))
custom_colors <- c("Bulk Soil" = "lightgreen", "Soil surrounding root" = "darkgreen")

# 
ggplot(abund_summary, aes(x = Family, y = mean_abund, fill = sample)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), color = "black") +
  geom_errorbar(aes(ymin = mean_abund - se_abund, ymax = mean_abund + se_abund),
                position = position_dodge(width = 0.9), width = 0.2, color = "black") +
  scale_fill_manual(values = custom_colors) +
  coord_flip() +
  labs(x = NULL, y = "Relative Abundance", fill = "Soil Compartment") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_text(face = "plain", color = "black"),
    axis.text.x = element_text(color = "black"),
    axis.title.y = element_text(color = "black"),
    axis.title.x = element_text(color = "black"),
    legend.text = element_text(color = "black"),
    legend.title = element_text(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.length = unit(-0.25, "cm"),
    axis.ticks = element_line(color = "black"),
    axis.text.x.top = element_text(margin = margin(t = 5)),
    axis.text.y.right = element_text(margin = margin(r = 5))
  )

#compo_UC_bact.jpg 600*700

#For MC
# ALDEx2_kw
t1 <- trans_diff$new(dataset = data_MC, method = "ALDEx2_kw", group = "sample", taxa_level = "Family")
t1$plot_diff_abund(use_number = 1:11, group_order = c("Bulk Soil","Soil surrounding root"))

data_MC$cal_abund()
abund_family <- data_MC$taxa_abund$Family
abund_long <- melt(as.matrix(abund_family))
colnames(abund_long) <- c("Family", "SampleID", "Abundance")
abund_long$Family <- sub(".*f__", "", abund_long$Family)

df_diff <- data.frame(
  Taxa = t1$res_diff$Taxa,
  P.adj = t1$res_diff$P.adj,
  stringsAsFactors = FALSE
)
df_diff <- df_diff[!grepl("f__NA", df_diff$Taxa), ]
df_diff$Family <- sub(".*f__", "", df_diff$Taxa)

top_families <- df_diff %>%
  arrange(P.adj) %>%
  slice(1:10) %>%
  pull(Family)

df_diff %>%
  filter(Family %in% top_families) %>%
  arrange(P.adj)
abund_long$Family <- factor(abund_long$Family, levels = rev(top_families))

abund_plot <- abund_long %>%
  filter(Family %in% top_families)

meta <- data_MC$sample_table
meta$SampleID <- rownames(meta)
meta
abund_plot <- abund_plot %>%
  left_join(meta, by = "SampleID")

abund_summary <- abund_plot %>%
  group_by(Family, sample) %>%
  summarise(
    mean_abund = mean(Abundance, na.rm = TRUE),
    se_abund = sd(Abundance, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

abund_summary$Family <- factor(abund_summary$Family, levels = levels(abund_plot$Family))
abund_summary$sample <- factor(abund_summary$sample, levels = c("Soil surrounding root","Bulk Soil"))
custom_colors <- c("Bulk Soil" = "yellow", "Soil surrounding root" = "gold")

# 
ggplot(abund_summary, aes(x = Family, y = mean_abund, fill = sample)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), color = "black") +
  geom_errorbar(aes(ymin = mean_abund - se_abund, ymax = mean_abund + se_abund),
                position = position_dodge(width = 0.9), width = 0.2, color = "black") +
  scale_fill_manual(values = custom_colors) +
  coord_flip() +
  labs(x = NULL, y = "Relative Abundance", fill = "Soil Compartment") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_text(face = "plain", color = "black"),
    axis.text.x = element_text(color = "black"),
    axis.title.y = element_text(color = "black"),
    axis.title.x = element_text(color = "black"),
    legend.text = element_text(color = "black"),
    legend.title = element_text(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.length = unit(-0.25, "cm"),
    axis.ticks = element_line(color = "black"),
    axis.text.x.top = element_text(margin = margin(t = 5)),
    axis.text.y.right = element_text(margin = margin(r = 5))
  )
#compo_MC_bact.jpg 600*700


#HC ALDEx2_kw
t1 <- trans_diff$new(dataset = data_HC, method = "ALDEx2_kw", group = "sample", taxa_level = "Family")
t1$plot_diff_abund(use_number = 1:11, group_order = c("Bulk Soil","Soil surrounding root"))

data_HC$cal_abund()
abund_family <- data_HC$taxa_abund$Family

abund_long <- melt(as.matrix(abund_family))
colnames(abund_long) <- c("Family", "SampleID", "Abundance")
abund_long$Family <- sub(".*f__", "", abund_long$Family)

df_diff <- data.frame(
  Taxa = t1$res_diff$Taxa,
  P.adj = t1$res_diff$P.adj,
  stringsAsFactors = FALSE
)

df_diff <- df_diff[!grepl("f__NA", df_diff$Taxa), ]
df_diff$Family <- sub(".*f__", "", df_diff$Taxa)

top_families <- df_diff %>%
  arrange(P.adj) %>%
  slice(1:10) %>%
  pull(Family)

df_diff %>%
  filter(Family %in% top_families) %>%
  arrange(P.adj)

abund_long$Family <- factor(abund_long$Family, levels = rev(top_families))

abund_plot <- abund_long %>%
  filter(Family %in% top_families)

meta <- data_HC$sample_table
meta$SampleID <- rownames(meta)
meta
abund_plot <- abund_plot %>%
  left_join(meta, by = "SampleID")

abund_summary <- abund_plot %>%
  group_by(Family, sample) %>%
  summarise(
    mean_abund = mean(Abundance, na.rm = TRUE),
    se_abund = sd(Abundance, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

abund_summary$Family <- factor(abund_summary$Family, levels = levels(abund_plot$Family))
abund_summary$sample <- factor(abund_summary$sample, levels = c("Soil surrounding root","Bulk Soil"))
custom_colors <- c("Bulk Soil" = "coral", "Soil surrounding root" = "red")

# 
ggplot(abund_summary, aes(x = Family, y = mean_abund, fill = sample)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), color = "black") +
  geom_errorbar(aes(ymin = mean_abund - se_abund, ymax = mean_abund + se_abund),
                position = position_dodge(width = 0.9), width = 0.2, color = "black") +
  scale_fill_manual(values = custom_colors) +
  coord_flip() +
  labs(x = NULL, y = "Relative Abundance", fill = "Soil Compartment") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_text(face = "plain", color = "black"),
    axis.text.x = element_text(color = "black"),
    axis.title.y = element_text(color = "black"),
    axis.title.x = element_text(color = "black"),
    legend.text = element_text(color = "black"),
    legend.title = element_text(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.length = unit(-0.25, "cm"),
    axis.ticks = element_line(color = "black"),
    axis.text.x.top = element_text(margin = margin(t = 5)),
    axis.text.y.right = element_text(margin = margin(r = 5))
  )
 }
###

# Functionality####
data_phylo_filt_rar
tax_data = data.frame(tax_table(data_phylo_filt_rar)) # create a separated file
otu_data = (data.frame(otu_table(data_phylo_filt_rar)))
sam_data = data.frame(sample_data(data_phylo_filt_rar))
otu_data_t <- as.data.frame(t(otu_data))
sam_data<-as.data.frame(sam_data)
rownames(sam_data)<-colnames(otu_data_t)
dataset <- microtable$new(sample_table = sam_data, otu_table = otu_data_t, tax_table = tax_data)
dataset

dataset
dataset$tax_table$Genus
dataset$tax_table$Specie <- dataset$tax_table$Genus
t1 <- trans_func$new(dataset)
t1$for_what<- 'prok'

t1$cal_spe_func(prok_database = "FAPROTAX" )

t1$res_spe_func
t1$cal_spe_func_perc(abundance_weighted = TRUE)
unique(colnames(t1$res_spe_func_perc))
df <- t1$res_spe_func_perc 
head(df)
head(dataset$sample_table)
dim(df)
dim(dataset$sample_table)
library(tidyverse)
sample_table <- dataset$sample_table %>%
  rownames_to_column("Sample")
if(!"Sample" %in% colnames(df)){
  df <- df %>%
    rownames_to_column("Sample")
}

df_joined <- df %>%
  left_join(sample_table, by = "Sample")
df_joined
df_long <- df_joined %>%
  pivot_longer(
    cols = c(chemoheterotrophy, nitrogen_fixation, ureolysis, denitrification),
    names_to = "Function",
    values_to = "Rel_Abundance"
  )

print(head(df_long))

df_avg <- df_long %>%
  group_by(Sample, sample, contamination, sampxcont, Function) %>%
  summarise(
    mean_perc = mean(Rel_Abundance, na.rm = TRUE),
    .groups = "drop"
  )

abund_summary <- df_avg %>%
  group_by(sample, contamination, sampxcont, Function) %>%
  summarise(
    mean_abund = mean(mean_perc, na.rm = TRUE),
    se_abund = sd(mean_perc, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )
print(abund_summary)

sample_colors_list <- list(
  "Uncont" = c("Bulk Soil" = "lightgreen", "Soil surrounding root" = "darkgreen"),
  "Mid Cont" = c("Bulk Soil" = "yellow", "Soil surrounding root" = "gold"),
  "High Cont" = c("Bulk Soil" = "coral", "Soil surrounding root" = "red")
)
cont_levels <- c("Uncont", "Mid Cont", "High Cont")

for (cont in cont_levels) {
  df_plot <- abund_summary %>%
    filter(contamination == cont) %>%
      mutate(sample = factor(sample, levels = c("Soil surrounding root","Bulk Soil")))
  custom_colors <- sample_colors_list[[cont]]
  library(forcats)

  p <- ggplot(df_plot, aes(x = fct_rev(Function), y = mean_abund, fill = sample)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9), color = "black") +
    geom_errorbar(aes(ymin = mean_abund - se_abund, ymax = mean_abund + se_abund),
                  position = position_dodge(width = 0.9), width = 0.2, color = "black") +
    scale_fill_manual(values = custom_colors) +
    coord_flip() +
    labs(
      title = paste(" "),
      x = NULL,
      y = "Relative Abundance (%)",
      fill = "Soil Compartment"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.y = element_text(face = "plain", color = "black"),
      axis.text.x = element_text(color = "black"),
      axis.title.y = element_text(color = "black"),
      axis.title.x = element_text(color = "black"),
      legend.text = element_text(color = "black"),
      legend.title = element_text(color = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks.length = unit(-0.25, "cm"),
      axis.ticks = element_line(color = "black"),
      axis.text.x.top = element_text(margin = margin(t = 5)),
      axis.text.y.right = element_text(margin = margin(r = 5)),
      plot.title = element_text(size = 14, face = "bold")
    )
  
  print(p)
}
#450*250
