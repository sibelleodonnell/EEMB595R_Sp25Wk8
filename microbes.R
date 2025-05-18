#Igraph walkthrough for network analysis
setwd('/Users/sibelleodonnell/Desktop/microbes')
install.packages("agricolae")
library(igraph)
library(igraphdata)
library(SpiecEasi)
library(microeco)
library(GUniFrac)
library(meconetcomp)
library(rgexf)
library(pheatmap)
library(aplot)
library(agricolae)
library(ggplot2)
library(magrittr)

# ------------- Microbial data wrangling with microtable (microeco)
#Working with built in 16S microbial soil data from An et al. 2019 
# surveyed soil prokaryotic communities in Chinese inland wetlands (IW), coastal wetland (CW) and Tibet plateau wetlands (TW) using 16S rRNA amplicon sequencing
# metadata table; data.frame
data(sample_info_16S)
# feature table; data.frame
data(otu_table_16S)
# taxonomic assignment table; data.frame
data(taxonomy_table_16S)
# phylogenetic tree
# Newick format; use read.tree function of ape package to read a tree
data(phylo_tree_16S)
# load the environmental data table if it is not in sample table
data(env_data_16S)

# fix the random number generation to make the results repeatable
set.seed(123)
# make the plotting background same with the tutorial
library(ggplot2)
theme_set(theme_bw())

#OTU table to dataframe
class(otu_table_16S)
otu_table_16S[1:5, 1:5]

#16S taxonomy table
#Do NOT use pure numbers as sample anmed in sample_table
class(taxonomy_table_16S)
taxonomy_table_16S[1:5, 1:3]
#If your data has NA's, unknowns, you can use tidy_taxonomy() function to clean
#tidy_taxonomy() also unifies taxonomic prefix automatically
taxonomy_table_16S %<>% tidy_taxonomy

#Sample metadata
class(sample_info_16S)
sample_info_16S[1:5, ]

#Environmental metadata stored separately (doesn't have to be)
class(env_data_16S)
env_data_16S[1:5, ]

#OTU phylogeny tree created in earlier analyses
class(phylo_tree_16S)

#Create an object of microtable class (similar to operation with package phyloseq)
class(mt)
# add data to microtable
mt <- microtable$new(sample_table = sample_info_16S, otu_table = otu_table_16S, tax_table = taxonomy_table_16S, phylo_tree = phylo_tree_16S)
mt


#Remove OTUs not assigned in kindgoms archaea or bacteria
mt$tax_table %<>% base::subset(Kingdom == "k__Archaea" | Kingdom == "k__Bacteria")
mt$tax_table %<>% .[grepl("Bacteria|Archaea", .$Kingdom), ]

#Remove OTUs assinged to mitochondria or chloroplasts
mt$filter_pollution(taxa = c("mitochondria", "chloroplast"))
mt

#Clean up data
mt$tidy_dataset()
mt

#To save data in microtable object to local files use save_table
mt$save_table(dirpath = "basic_files", sep = ",")

#Calculate taxonomic abundance at each rank using cal_abund()
#Can use this function for both relative and absolute abundance calculations
mt$cal_abund()
class(mt$taxa_abund)
mt$taxa_abund$Phylum[1:5, 1:5]

#Resampling to reduce impact of sequencing depth on diversity measurements
mt_rarefied <- clone(mt)
# use sample_sums to check the sequence numbers in each sample
mt_rarefied$sample_sums() %>% range

#Using 10000 sequences per sample
mt_rarefied$rarefy_samples(sample.size = 10000)
mt_rarefied$sample_sums() %>% range

#Calculate alpha diversity
mt_rarefied$cal_alphadiv(PD = FALSE)
#Return alpha diversity in object
class(mt_rarefied$alpha_diversity)
# save alpha_diversity
mt_rarefied$save_alphadiv(dirpath = "alpha_diversity")

#Calulate beta diversity
#If method parameter is not provided, the function automatically calculates Bray-curtis, Jaccard, weighted Unifrac and unweighted unifrac matrixes
mt_rarefied$cal_betadiv(unifrac = TRUE)
# return beta_diversity list in the object
class(mt_rarefied$beta_diversity)
# save beta_diversity to a directory
mt_rarefied$save_betadiv(dirpath = "beta_diversity")

#With microtable you can also merge, subset by sample or taxa, rename features in all files of microtable object, etc etc etc
#TLDR: microtable rocks!

# ------------- Mic
data(soil_amp)

#Create a list
soil_amp_network <- list()
#Save raw data
tmp <- clone(soil_amp)

#Select samples from IW group and change sample_table directly
tmp$sample_table %<>% subset(Group == "IW")
#trimming with tidy_dataset()
tmp$tidy_dataset()
#filtering out low relative abundance
tmp <- trans_network$new(dataset = tmp, cor_method = "spearman", filter_thres = 0.0005)
#p-value and correlation coefficient thresholds
tmp$cal_network(COR_p_thres = 0.01, COR_cut = 0.6)

#Network put into a list 
soil_amp_network$IW <- tmp

#Repeat steps above for TW and CW groups
# select samples of "TW" group
tmp <- clone(soil_amp)
tmp$sample_table %<>% subset(Group == "TW")
tmp$tidy_dataset()
tmp <- trans_network$new(dataset = tmp, cor_method = "spearman", filter_thres = 0.0005)
tmp$cal_network(COR_p_thres = 0.01, COR_cut = 0.6)
soil_amp_network$TW <- tmp
# select samples of "CW" group
tmp <- clone(soil_amp)
tmp$sample_table %<>% subset(Group == "CW")
tmp$tidy_dataset()
tmp <- trans_network$new(dataset = tmp, cor_method = "spearman", filter_thres = 0.0005)
tmp$cal_network(COR_p_thres = 0.01, COR_cut = 0.6)
soil_amp_network$CW <- tmp
#Yay our network is now a list by each soil type

#Partitioning modules for all networks in the list (IW,CW,TW)
soil_amp_network %<>% cal_module(undirected_method = "cluster_fast_greedy")
# "cluster_fast_greedy" is an algorithm from the igraph package that detects communities in undirected graphs by greedily optimizing modularity
  # Undirected means that the connections (edges) dont have a direction
# a module of microbes is a community more connected to eachother than those outside of the module 
#Save topological attributes for the networks (please don't ask me to explain what these mean...)
tmp <- cal_network_attr(soil_amp_network)
#Extract node and edge properties of networks
soil_amp_network %<>% get_node_table(node_roles = TRUE) %>% get_edge_table

#Compare nodes shared across networks
tmp <- node_comp(soil_amp_network, property = "name")
# obtain nodes intersection
tmp1 <- trans_venn$new(tmp, ratio = "numratio")
g1 <- tmp1$plot_venn(fill_color = FALSE)
quartz()
print(g1)
# calculate jaccard distance to reflect the overall differences of networks
tmp$cal_betadiv(method = "jaccard")
tmp$beta_diversity$jaccard

#Compare edges across networks
tmp <- edge_comp(soil_amp_network)
# obtain edges intersection
tmp1 <- trans_venn$new(tmp, ratio = "numratio")
g2 <- tmp1$plot_venn(fill_color = FALSE)
quartz()
print(g2)
# calculate jaccard distance
tmp$cal_betadiv(method = "jaccard")
tmp$beta_diversity$jaccard

#Extract overlapped edges of networks to a new network
# first obtain edges distribution and intersection
tmp <- edge_comp(soil_amp_network)
tmp1 <- trans_venn$new(tmp)
# convert intersection result to a microtable object
tmp2 <- tmp1$trans_comm()
# extract the intersection of all the three networks ("IW", "TW" and "CW")
#use colnames(tmp2$otu_table) to find the required name
Intersec_all <- subset_network(soil_amp_network, venn = tmp2, name = "IW&TW&CW")
# Intersec_all is a trans_network object
# for example, save Intersec_all as gexf format
Intersec_all$save_network("Intersec_all.gexf")



###### Walkthrough plotting
#### This code below is essentially doing all that we did above in one shot ####
test <- microtable$
  new(sample_table = sample_info_16S, otu_table = otu_table_16S, tax_table = taxonomy_table_16S, phylo_tree = phylo_tree_16S)$
  filter_pollution(taxa = c("mitochondria", "chloroplast"))$
  tidy_dataset()$
  filter_taxa()

#Testing sequencing depth per sample
test$sample_sums() %>% range
#min of 10341, max of 10341 reads per sample (practice data)

#### Same thing, we did this earlier ####
test$
  rarefy_samples()$ #subsampling to an even depth for all samples for fair diversity comparison
  cal_abund()$ #calculate relative abundance
  save_abund(file.path("taxa_abund"))$
  cal_alphadiv()$ #alpha diversity calculation
  save_alphadiv(file.path("alpha_diversity"))$
  cal_betadiv()$ #beta diversity calculations
  save_betadiv(file.path("beta_diversity"))



# trans_beta

library(magrittr)
library(microeco)
data(dataset)

#Computing PCoA ordination basd on Bray-Curtis distances (why use bray curtis distances?)
t1 <- trans_beta$
  new(dataset = dataset, group = "Group", measure = "bray")$
  cal_ordination(method = "PCoA", ncomp = 5)

g1 <- t1$plot_ordination(plot_color = "Group")
quartz()
g1
#^Each point is a sample from one of the three regions. Closer points have more similar microbial structures

trans_beta$
  new(dataset = dataset, group = "Type", measure = "bray")$
  cal_group_distance(within_group = TRUE)$
  cal_group_distance_diff(method = "wilcox")$
  plot_group_distance()
#^ within-group vs. between-group dissimilarities (beta diversity) 
  #IW: NE, NW 
  #CW: NC, YML, SC
  #TW: QTP 

trans_beta$
  new(dataset = dataset, measure = "bray")$
  cal_manova(manova_all = FALSE, group = "Type", by_group = "Group")$
  res_manova
#Some significant differences between types (sites) within a group (region)

trans_beta$
  new(dataset = dataset, group = "Group", measure = "bray")$
  cal_anosim()$
  res_anosim
#Analysis of similarities non-parametric (does not assume data follows a specific distribution; used for not normally distributed data) test for group separation 


# Microbial network analysis (Finally)

library(microeco)
data(dataset)

#This plot is so incredibly crowded and messy (fix it)
t1<- trans_network$new(dataset = dataset, cor_method = "pearson", filter_thres = 0.001)$
  cal_network(COR_p_thres = 0.05, COR_cut = 0.6)$ #filtering thresholds for which edges get retained
  cal_module()$
  cal_network_attr()$
  get_node_table()$
  get_edge_table()$
  cal_sum_links()$
  save_network()
quartz()
t1$plot_network()
#^Ew
#Pretty circular chord diagram
quartz()
t1$plot_sum_links(method = "circlize", transparency = 0.2, annotationTrackHeight = circlize::mm_h(c(5, 5)))
# Groups with wider arcs are more interaction-dense — more co-occurrence.
#Chord coloring is sneaky here and doesn't mean it has any directionality (we cannot determine that from microbial community analyses like these -- need time series or intervention data)
  #In this case, color just starts at the black node as priority color and then loops around counter clockwise
# Chords that loop back to the same arch mean that there are taxa that have  strong co-occurances and might be co-adapted
quartz()
t1$plot_taxa_roles(use_type = 1)
















