#################################################################################
#
#           Gut microbiome structure and function in wild spotted hyenas
#                      
#           Rojas et al 2022. Taxonomic and functional variation of the gut 
#               microbiome over host lifespan in wild spotted hyenas
#
#                         Created By: Cassie Ettinger, & Jason Stajich. (2022)
#                         Adapted By: Connie Rojas
#                         Last updated:  25 Jul 2022
#
################################################################################

## CODE FOR: making phylogenetic trees of individual MAGs, placed within the 
#             broader GTDB-r95 taxonomy

source(file="scripts/00_background.R"); #load necessary packages and specifications


################################################################################
#         1.  Load taxonomic and phylogenetic data
################################################################################

# taxonomic classifications of your bins as determined by GTDB release 95
bin_names=read.csv("data/metagenomes/00_GTDB_MAG_taxonomy_for_trees.csv", 
                   header=T);

# giant phylogeny output by GTDB when assigning taxonomy to MAGs
phylo.r95 <- read.newick("data/metagenomes/00_gtdbtk.bac120.classify.tree");
phylo.r95$edge.length <- phylo.r95$edge.length*10;

# all GTDB-r95 taxonomic labels downloaded from https://data.gtdb.ecogenomic.org/releases/
load("data/metagenomes/00_GTDB_r95_taxonomy.Rdata");


################################################################################
#         2.  Subset GTDB phylogeny to only MAG of interest and replace 
#                 tree node labels with GTDB taxonomic labels
################################################################################

# select MAG of interest and how how many node levels to show
# 4 - broad while 2 - narrow; 
mag.tree <- tree_subset(phylo.r95, "metabat.306", levels_back = 4);

# add your MAG taxonomic labels to the larger GTDB taxonomy
tree.meta.ms <- full_join(taxonomy.bac, bin_names);

# subset larger GTDB taxonomy to only your labels of interest
mag.tree.taxaIDs <- mag.tree$tip.label;
tree.meta.ms2 <- tree.meta.ms %>% filter(label %in% mag.tree.taxaIDs);

# add taxonomic labels to your tree in preparation for plotting
my.tree <-full_join(mag.tree, tree.meta.ms2, by = "label");


################################################################################
#     3.  Plot phylogeny of MAG placed within the greater GTDB phylogeny
################################################################################

# plot tree
ptree=ggtree(my.tree, color = "black", size = 1, linetype = 1) + 
  geom_tiplab(aes(label = Taxonomy, color=Source),
              fontface = "bold.italic", size = 2.5, offset = 0.001) +
  xlim(0,6) + # (0,11) if detailed tree
  scale_color_manual(values =c(study = "#EF7F4FFF", ref = "#000000")) +
  theme(legend.position = "none"); plot(ptree);

# save tree
ggsave(filename = 'figures/metagenomes/07_MAG306_phylogeny.pdf', 
       plot = last_plot(), 
       device = 'pdf', 
       width = 10, 
       height = 7, 
       dpi = 300);














