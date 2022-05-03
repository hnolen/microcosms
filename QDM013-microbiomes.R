#####QDM-013 microbiome work

#####ITS work
'''
everything so far has been done in QIIME2-2020.2
'''

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("phyloseq")

install.packages("remotes")
remotes::install_github("jbisanz/qiime2R")

library(phyloseq)
library(qiime2R)
library(ggplot2)
library(ape)
library(vegan)
library(ggpubr)
library(agricolae)
library(ggsignif)
library(dplyr)


##following code from the QIIME2_to_phyloseq.pdf
#i didn't filter anything out but doing these steps just to make sure everything is good for the phyloseq step
#need to merge these files in R and output a merged file
#read in feature table
otu<-read.table(file = "ITS/otu_table.txt", header = TRUE, check.names = FALSE)
head(otu)

tax<-read.table(file = "ITS/taxonomy.tsv", header = TRUE, sep = "\t")
head(tax)

#merge files
merged_file<-merge(otu, tax, by.x = c("OTUID"), by.y = c("OTUID"))
head(merged_file)

#note: number of rows should equal your shortest file length, drops taxonomy for OTUs that don't exist in your feature table

#output merged .txt file
write.table(merged_file, file = "./ITS/combined_otu_tax.txt", sep = "\t", col.names = TRUE, row.names = FALSE)

'''
Need to open the merged .txt file in excel and split into two files: one for taxonomy
(with columns OTUID and taxonmic info) and the other file for the feature table (OTUID and abundances)
save each of these as a .csv

saved as split_feat_table.csv and split-its-tax.csv
'''


###### ITS data

feat_table = read.csv("ITS/split_feat_table.csv", sep = ",", row.names =1, check.names = FALSE)
feat_table = as.matrix(feat_table)

taxonomy = read.csv("ITS/split-its-tax.csv", sep = ",", row.names = 1)
taxonomy = as.matrix(taxonomy)

metadata = read.table("ITS/metadata.tsv", sep = "\t", row.names =1, header = TRUE)

phy_tree = read_tree("ITS/tree.nwk")

#import as phyloseq objects
OTU = otu_table(feat_table, taxa_are_rows = TRUE)
TAX = tax_table(taxonomy)
META = sample_data(metadata)
#tree was already imported as phyloseq object


#https://joey711.github.io/phyloseq/import-data.html

physeq = phyloseq(OTU, TAX, META, phy_tree)

##bar plots - https://joey711.github.io/phyloseq/plot_bar-examples.html 
plot_bar(physeq, x = "host_inoculated", y = "Abundance", fill = "Phylum")

##heat maps - https://joey711.github.io/phyloseq/plot_heatmap-examples.html 
plot_heatmap(physeq,sample.label = "microbiome_origin", taxa.label="Phylum")

##trying to get ordination plot - https://joey711.github.io/phyloseq/plot_ordination-examples.html#mds_(%E2%80%9Cpcoa%E2%80%9D)_on_unifrac_distances
ordu = ordinate(physeq, "PCoA", "bray") #other options like "NMDS" and "unifrac"
p1<-plot_ordination(physeq, ordu, color = "host_inoculated", shape = "microbiome_origin") +
  geom_point(size = 2) +
  theme(axis.title.x = element_text(color="black", size = 15)) +
  theme(axis.title.y = element_text(color="black", size =15))

#plot microbiome network
plot_net(physeq, maxdist = 0.4, point_label = "host_inoculated", color = "host_inoculated", shape = "microbiome_origin")

'''
https://joey711.github.io/phyloseq/plot_network-examples.html 
can change the microbiome network by doing:
ig<-make_network(physeq1, dist.fun = "bray", max.dist = 0.3)
and then using the ig as the data object in the plot_net() function - plot_net defaults to jaccard distance
'''


#richness plots - https://joey711.github.io/phyloseq/plot_richness-examples.html 
p<-plot_richness(physeq, x = "host_inoculated", measures = c("Chao1", "Shannon"))
p +
  geom_boxplot()

p<-plot_richness(physeq, x = "microbiome_origin", measures = c("Chao1", "Shannon"))
p +
  geom_boxplot()
