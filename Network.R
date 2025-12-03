#Scripts for the article "Microbial refuge under shrubs limits soil Cu contamination effects"
#Submitted to Soil Biology and Biochemistry
#
#December 2nd, 2025

# Inquiries to Dr. Pedro Mondaca > pedro.mondaca@pucv.cl


{
library(dplyr)
library(igraph)
library(ggraph)
library(tidygraph)
library(WGCNA)
library(tibble)
library(reshape2)
library(microeco) # for optional functional annotation
library(NetCoMi)
library(phyloseq)
library(network)
library(sna)
library(tidyverse)
library(tidygraph)
library(ggClusterNet)
library(datastorr)
library(fungaltraits)
library(scales)
library(ggrepel)
library(graphlayouts)

data_phylo_filt_rar_fam_b
data_phylo_filt_rar_fam_f

phylo_bac_f
phylo_fun_f

tax_table(phylo_bac_f)
tax_table(phylo_fun_f)

tax_bf<-tax_table(phylo_bac_f)
tax_ff<-tax_table(phylo_fun_f)
tax_bf<-tax_bf[,-c(6:7)]
tax_ff<-tax_ff[,-c(6:7)]
tax_bf[,1:5]
tax_ff[,1:5]

otu_bf<-otu_table(phylo_bac_f)
otu_ff<-otu_table(phylo_fun_f)
otu_bf[,1:3]
otu_ff[,1:3]

sam_bf<-sample_data(phylo_bac_f)
sam_ff<-sample_data(phylo_fun_f)
sam_bf
sam_ff

otu_ff_reordered <- otu_ff[match(rownames(otu_bf), rownames(otu_ff)), ]
sam_ff_reordered <- sam_ff[match(rownames(sam_bf), rownames(sam_ff)), ]

all(rownames(otu_ff_reordered) == rownames(otu_bf))  # should be TRUE

OTU = otu_table(as.matrix(otu_ff_reordered), taxa_are_rows = FALSE) # create the occurrence table object in phyloseq format
SAM = sample_data(sam_ff_reordered, errorIfNULL = TRUE) # create the sample metadata object in phyloseq format
TAX = tax_table(as.matrix(tax_ff)) # create the observation metadata object (OTU taxonomy) in phyloseq format
phylo_ff_reordered <- phyloseq(OTU, TAX, SAM) # create the phyloseq object including occurrence table data and sample/observation metadata
phylo_ff_reordered # print information about the phyloseq object


OTU = otu_table(as.matrix(otu_bf), taxa_are_rows = FALSE) # create the occurrence table object in phyloseq format
SAM = sample_data(sam_bf, errorIfNULL = TRUE) # create the sample metadata object in phyloseq format
TAX = tax_table(as.matrix(tax_bf)) # create the observation metadata object (OTU taxonomy) in phyloseq format
phylo_bf_reordered <- phyloseq(OTU, TAX, SAM) # create the phyloseq object including occurrence table data and sample/observation metadata
phylo_bf_reordered # print information about the phyloseq object

#unamos bact y fungi cata
ps16 = phylo_bf_reordered
psIT = phylo_ff_reordered
ps16
psIT
taxa_names(ps16)
taxa_names(psIT)
tax_table(ps16)
tax_table(psIT)

colnames(tax_table(ps16)) <- c("Kingdom", "Phylum", "Class", "Order", "Family")
colnames(tax_table(psIT)) <- c("Kingdom", "Phylum", "Class", "Order", "Family")


sample_names(ps16)
sample_names(psIT)
all(sample_names(ps16) == sample_names(psIT))
colnames(sample_data(ps16))
colnames(sample_data(psIT))
rownames(sample_data(ps16))
rownames(sample_data(psIT))
sapply(sample_data(ps16), class)
sapply(sample_data(psIT), class)
sample_data(psIT) <- sample_data(ps16)
rownames(sample_data(psIT)) == rownames(sample_data(ps16))

ps = ggClusterNet::merge16S_ITS(ps16s = ps16,
                                psITS = psIT,
                                N16s = 136,
                                NITS = 145,
                                scale = TRUE,
                                onlygroup = FALSE,
                                dat1.lab = "bac",
                                dat2.lab = "fun"
)
ps
sample_data(ps)
tax_table(ps)


# Subset by sampxcont
ps_BS_UC <- subset_samples(ps, sampxcont == "BS UC")
ps_BS_UC <- prune_taxa(taxa_sums(ps_BS_UC) > 0, ps_BS_UC)
otu_BS_UC <- as.data.frame(otu_table(ps_BS_UC))

ps_BS_MC <- subset_samples(ps, sampxcont == "BS MC")
ps_BS_MC <- prune_taxa(taxa_sums(ps_BS_MC) > 0, ps_BS_MC)
otu_BS_MC <- as.data.frame(otu_table(ps_BS_MC))

ps_BS_HC <- subset_samples(ps, sampxcont == "BS HC")
ps_BS_HC <- prune_taxa(taxa_sums(ps_BS_HC) > 0, ps_BS_HC)
otu_BS_HC <- as.data.frame(otu_table(ps_BS_HC))

ps_SSR_UC <- subset_samples(ps, sampxcont == "SSR UC")
ps_SSR_UC <- prune_taxa(taxa_sums(ps_SSR_UC) > 0, ps_SSR_UC)
otu_SSR_UC <- as.data.frame(otu_table(ps_SSR_UC))

ps_SSR_MC <- subset_samples(ps, sampxcont == "SSR MC")
ps_SSR_MC <- prune_taxa(taxa_sums(ps_SSR_MC) > 0, ps_SSR_MC)
otu_SSR_MC <- as.data.frame(otu_table(ps_SSR_MC))

ps_SSR_HC <- subset_samples(ps, sampxcont == "SSR HC")
ps_SSR_HC <- prune_taxa(taxa_sums(ps_SSR_HC) > 0, ps_SSR_HC)
otu_SSR_HC <- as.data.frame(otu_table(ps_SSR_HC))

dim(otu_BS_UC)
otu_BS_UC2 <-otu_BS_UC + 1e-6
otu_BS_UC2
colSums(otu_BS_UC2)
otu_BS_UC3 <-  sweep(otu_BS_UC, 2, colSums(otu_BS_UC), FUN = "/")
colSums(otu_BS_UC3)
otu_BS_UC3
otu_BS_UC4 <- t(otu_BS_UC3)
dim(otu_BS_UC4)
otu_BS_UC4
rowSums(otu_BS_UC4)

dim(otu_BS_MC)
otu_BS_MC2 <-otu_BS_MC + 1e-6
otu_BS_MC2
colSums(otu_BS_MC2)
otu_BS_MC3 <-  sweep(otu_BS_MC, 2, colSums(otu_BS_MC), FUN = "/")
colSums(otu_BS_MC3)
otu_BS_MC3
otu_BS_MC4 <- t(otu_BS_MC3)
dim(otu_BS_MC4)
otu_BS_MC4
rowSums(otu_BS_MC4)

dim(otu_BS_HC)
otu_BS_HC2 <-otu_BS_HC + 1e-6
otu_BS_HC2
colSums(otu_BS_HC2)
otu_BS_HC3 <-  sweep(otu_BS_HC, 2, colSums(otu_BS_HC), FUN = "/")
colSums(otu_BS_HC3)
otu_BS_HC3
otu_BS_HC4 <- t(otu_BS_HC3)
dim(otu_BS_HC4)
otu_BS_HC4
rowSums(otu_BS_HC4)

dim(otu_SSR_UC)
otu_SSR_UC2 <-otu_SSR_UC + 1e-6
otu_SSR_UC2
colSums(otu_SSR_UC2)
otu_SSR_UC3 <-  sweep(otu_SSR_UC, 2, colSums(otu_SSR_UC), FUN = "/")
colSums(otu_SSR_UC3)
otu_SSR_UC3
otu_SSR_UC4 <- t(otu_SSR_UC3)
dim(otu_SSR_UC4)
otu_SSR_UC4
rowSums(otu_SSR_UC4)

dim(otu_SSR_MC)
otu_SSR_MC2 <-otu_SSR_MC + 1e-6
otu_SSR_MC2
colSums(otu_SSR_MC2)
otu_SSR_MC3 <-  sweep(otu_SSR_MC, 2, colSums(otu_SSR_MC), FUN = "/")
colSums(otu_SSR_MC3)
otu_SSR_MC3
otu_SSR_MC4 <- t(otu_SSR_MC3)
dim(otu_SSR_MC4)
otu_SSR_MC4
rowSums(otu_SSR_MC4)

dim(otu_SSR_HC)
otu_SSR_HC2 <-otu_SSR_HC + 1e-6
otu_SSR_HC2
colSums(otu_SSR_HC2)
otu_SSR_HC3 <-  sweep(otu_SSR_HC, 2, colSums(otu_SSR_HC), FUN = "/")
colSums(otu_SSR_HC3)
otu_SSR_HC3
otu_SSR_HC4 <- t(otu_SSR_HC3)
dim(otu_SSR_HC4)
otu_SSR_HC4
rowSums(otu_SSR_HC4)

}


#BS UC
{
fVar <- 250 #n° of nodes
net_BS_UC <- netConstruct(
  otu_BS_UC4,
  measure     = "spieceasi",
  normMethod  = "none",
  zeroMethod  = "pseudo",    
  filtTax     = "highestVar",
  filtTaxPar  = list(highestVar = fVar),
  
  # SpiecEasi parameters
  measurePar  = list(
    method           = "mb",        
    sel.criterion    = "stars",    
    pulsar.params    = list(
      rep.num          = 30,       
      thresh           = 0.05,     
      subsample.ratio  = 0.8,
      seed             = 2015
    ),
    nlambda          = 20,         
    lambda.min.ratio = 0.05        
  ),
  sparsMethod = "none",
  verbose     = 3,
  seed        = 2015
)

net_BS_UC_anal <- netAnalyze(net_BS_UC)
summary(net_BS_UC_anal)

edgelist <- net_BS_UC$edgelist1[
  order(net_BS_UC$edgelist1$adja, decreasing = TRUE), ]
edgelist |> head()
edgelist |> tail()

netcomi_graph <- SpiecEasi::adj2igraph(abs(net_BS_UC$adjaMat1))
netcomi_graph

str(net_BS_UC_anal)

# Get edgelist 
edg <- net_BS_UC$edgelist1 %>%
  mutate(sign = ifelse(asso >= 0, "positive", "negative"),
         w_abs = abs(asso))

g <- graph_from_data_frame(edg, directed = FALSE)

V(g)$type <- ifelse(grepl("^bac_", V(g)$name), "Bacteria",
                    ifelse(grepl("^fun_", V(g)$name), "Fungi", "Unknown"))
V(g)$deg  <- igraph::degree(g) #degree
V(g)$eig  <- igraph::eigen_centrality(g)$vector #eigenvector centrality
V(g)$comm <- cluster_fast_greedy(g)$membership   # or cluster_louvain(g)

# Edge aesthetics
E(g)$weight <- edg$w_abs
E(g)$sign   <- edg$sign
E(g)$w_abs <- ifelse(is.na(E(g)$weight), 0, abs(E(g)$weight))

hub_n  <- 10
hub_ids <- order(V(g)$eig, decreasing = TRUE)[seq_len(min(hub_n, gorder(g)))]
V(g)$hub <- ifelse(seq_along(V(g)) %in% hub_ids, V(g)$name, NA)

if (is.null(V(g)$eig)) V(g)$eig <- igraph::evcent(g, scale = TRUE)$vector
if (is.null(V(g)$type)) {
  V(g)$type <- ifelse(grepl("^bac_", V(g)$name), "Bacteria",
                      ifelse(grepl("^fun_", V(g)$name), "Fungi", "Unknown"))
}
if (is.null(E(g)$sign)) {
  if (!is.null(E(g)$weight)) {
    E(g)$sign <- ifelse(E(g)$weight >= 0, "positive", "negative")
  } else stop("Falta 'sign' o 'weight' en edges para clasificar asociaciones.")
}

table(V(g)$type)


# ---------- Top 10 associations ----------
ed_df <- igraph::as_data_frame(g, what = "edges") %>%
  transmute(from, to, sign = as.character(sign))

count_by_sign <- ed_df %>%
  pivot_longer(c(from, to), names_to = "end", values_to = "node") %>%
  count(node, sign) %>%
  pivot_wider(names_from = sign, values_from = n, values_fill = 0) %>%
  rename(neg = negative, pos = positive)

top_pos_nodes <- count_by_sign %>% arrange(desc(pos)) %>% slice_head(n = 5) %>% pull(node)
top_neg_nodes <- count_by_sign %>% arrange(desc(neg)) %>% slice_head(n = 5) %>% pull(node)
top_pos_nodes
V(g)$label <- ifelse(V(g)$name %in% union(top_pos_nodes, top_neg_nodes), V(g)$name, NA)
V(g)$label_sign <- case_when(
  V(g)$name %in% top_pos_nodes & V(g)$name %in% top_neg_nodes ~ "both",
  V(g)$name %in% top_pos_nodes ~ "positive",
  V(g)$name %in% top_neg_nodes ~ "negative",
  TRUE ~ NA_character_
)


# Filter weak edges
g_filt <- igraph::delete_edges(
  g,
  E(g)[abs(weight) < 0.05]
)
V(g_filt)$degree <- igraph::degree(g_filt)

V(g_filt)$size <- scales::rescale(V(g_filt)$degree, to = c(1, 12))

set.seed(42)

xy <- layout_components(g_filt, layout = layout_with_stress)
lay <- create_layout(g_filt, layout = xy)
lay_df <- as.data.frame(lay) %>% filter(!is.na(label))

V(g_filt)$deg <- igraph::degree(g_filt, mode = "all", loops = FALSE)
node_size <- scales::rescale(V(g_filt)$deg, to = c(3, 12))   # <--- AQUÍ se define

p <- ggraph(lay) +
  geom_edge_link(aes(color = sign,
                     edge_width = abs(weight),
                     edge_alpha = abs(weight)),
                 show.legend = TRUE) +
  scale_edge_color_manual(values = c(positive = "#1f78b4", negative = "#e31a1c"),
                          name = "Asociación") +
  scale_edge_width(range = c(0.2, 1.8)) +     
  scale_edge_alpha(range = c(0.5, 0.8)) +    
  guides(edge_width = "none", edge_alpha = "none") +
  
  geom_node_point(aes(fill = type),
                  shape = 21, size = node_size, color = "grey25", stroke = 0.4) +
  scale_fill_manual(values = c(Bacteria = "#E69F00", Fungi = "#0072B2", Unknown = "grey70"),
                    name = "Dominio") +
  theme_void() +
  coord_equal()

p #900*900


if (is.null(E(g)$sign)) {
  if (is.null(E(g)$weight)) stop("La red no tiene 'sign' ni 'weight' en las aristas.")
  E(g)$sign <- ifelse(E(g)$weight >= 0, "positive", "negative")
}

if (is.null(V(g)$type)) {
  V(g)$type <- ifelse(grepl("^bac_", V(g)$name), "Bacteria",
                      ifelse(grepl("^fun_", V(g)$name), "Fungi", "Unknown"))
}

edge_df <- igraph::as_data_frame(g, what = "edges") %>%
  mutate(
    from_type = V(g)$type[match(from, V(g)$name)],
    to_type   = V(g)$type[match(to,   V(g)$name)],
    # interaction type
    interaction = dplyr::case_when(
      from_type == "Bacteria" & to_type == "Bacteria" ~ "B-B",
      from_type == "Fungi"    & to_type == "Fungi"    ~ "F-F",
      TRUE                                            ~ "B-F"
    ),
    association = ifelse(as.character(sign) %in% c("positive","Positive"), "Positive", "Negative")
  )

prop_df <- edge_df %>%
  group_by(interaction, association) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(interaction) %>%
  mutate(proportion = count / sum(count))

edge_counts <- edge_df %>%
  count(interaction, association, name = "n") %>%
  group_by(interaction) %>%
  mutate(
    percentage = n / sum(n) * 100,
    label = paste0(n, "\n(", round(percentage, 1), "%)")
  ) %>%
  ungroup()

total_counts <- edge_counts %>%
  group_by(interaction) %>%
  summarise(total = sum(n), .groups = "drop")

edge_counts <- left_join(edge_counts, total_counts, by = "interaction")
assoc_colors <- c("Positive" = "#56B4E9",  
                  "Negative" = "#D55E00")  

#plot associations
y_max <- max(edge_counts$total, na.rm = TRUE)
neg_offset <- 0.05 * y_max  # 

ggplot(edge_counts, aes(x = interaction, y = n, fill = association)) +
  geom_bar(stat = "identity", color = "black") +
  

  geom_text(
    data = dplyr::filter(edge_counts, association == "Positive"),
    aes(label = label),
    position = position_stack(vjust = 0.5),
    color = "black", size = 4
  ) +

  geom_text(
    data = dplyr::filter(edge_counts, association == "Negative"),
    aes(y = total + neg_offset, label = label),
    color = "black", size = 4
  ) +
  
  scale_fill_manual(values = assoc_colors) +
  labs(x = "", y = "Number of Edges", fill = "Association") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    axis.line        = element_line(color = "black")
  )
#493*536

ggplot(edge_counts, aes(x = interaction, y = n, fill = association)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = assoc_colors) +
  labs(x = "", y = "Number of Edges", fill = "Association") +
  scale_y_continuous(limits = c(0, 200)) +   # <-- this controls Y-axis limit
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    axis.line        = element_line(color = "black")
  )
#493*536

#Obtain negative edges
negative_ff <- edge_df %>%
  filter(interaction == "F-F", association == "Negative")
negative_ff

negative_bf <- edge_df %>%
  filter(interaction == "B-F", association == "Negative")
negative_bf

positive_bf <- edge_df %>%
  filter(interaction == "B-F", association == "Positive")
positive_bf

negative_bb <- edge_df %>%
  filter(interaction == "B-B", association == "Negative")
negative_bb


}


#BS MC
{ 
net_BS_MC <- netConstruct(
  otu_BS_MC4,
  measure     = "spieceasi",
  normMethod  = "none",
  zeroMethod  = "pseudo",    
  filtTax     = "highestVar",
  filtTaxPar  = list(highestVar = fVar),
   
   #  SpiecEasi parameters
  measurePar  = list(
    method           = "mb",        
    sel.criterion    = "stars",   
    pulsar.params    = list(
      rep.num          = 30,        
      thresh           = 0.05,      
      subsample.ratio  = 0.8,
      seed             = 2015
    ),
    nlambda          = 20,          
    lambda.min.ratio = 0.05         
  ),
   verbose     = 3,
  seed        = 2015
)

net_BS_MC_anal <- netAnalyze(net_BS_MC)
summary(net_BS_MC_anal)

edgelist <- net_BS_MC$edgelist1[
  order(net_BS_MC$edgelist1$adja, decreasing = TRUE), ]
edgelist |> head()
edgelist |> tail()

netcomi_graph <- SpiecEasi::adj2igraph(abs(net_BS_MC$adjaMat1))
netcomi_graph

str(net_BS_MC_anal)

edg <- net_BS_MC$edgelist1 %>%
  mutate(sign = ifelse(asso >= 0, "positive", "negative"),
         w_abs = abs(asso))

g <- graph_from_data_frame(edg, directed = FALSE)

V(g)$type <- ifelse(grepl("^bac_", V(g)$name), "Bacteria",
                    ifelse(grepl("^fun_", V(g)$name), "Fungi", "Unknown"))
V(g)$deg  <- igraph::degree(g) #degree
V(g)$eig  <- igraph::eigen_centrality(g)$vector #eigenvector centrality
V(g)$comm <- cluster_fast_greedy(g)$membership   # or cluster_louvain(g)

# Edge aesthetics
E(g)$weight <- edg$w_abs
E(g)$sign   <- edg$sign
E(g)$w_abs <- ifelse(is.na(E(g)$weight), 0, abs(E(g)$weight))

hub_n  <- 10
hub_ids <- order(V(g)$eig, decreasing = TRUE)[seq_len(min(hub_n, gorder(g)))]
V(g)$hub <- ifelse(seq_along(V(g)) %in% hub_ids, V(g)$name, NA)

if (is.null(V(g)$eig)) V(g)$eig <- igraph::evcent(g, scale = TRUE)$vector
if (is.null(V(g)$type)) {
  V(g)$type <- ifelse(grepl("^bac_", V(g)$name), "Bacteria",
                      ifelse(grepl("^fun_", V(g)$name), "Fungi", "Unknown"))
}
if (is.null(E(g)$sign)) {
  if (!is.null(E(g)$weight)) {
    E(g)$sign <- ifelse(E(g)$weight >= 0, "positive", "negative")
  } else stop("Falta 'sign' o 'weight' en edges para clasificar asociaciones.")
}

table(V(g)$type)


# ---------- Top 10 associations ----------
ed_df <- igraph::as_data_frame(g, what = "edges") %>%
  transmute(from, to, sign = as.character(sign))

count_by_sign <- ed_df %>%
  pivot_longer(c(from, to), names_to = "end", values_to = "node") %>%
  count(node, sign) %>%
  pivot_wider(names_from = sign, values_from = n, values_fill = 0) %>%
  rename(neg = negative, pos = positive)

top_pos_nodes <- count_by_sign %>% arrange(desc(pos)) %>% slice_head(n = 5) %>% pull(node)
top_neg_nodes <- count_by_sign %>% arrange(desc(neg)) %>% slice_head(n = 5) %>% pull(node)
top_pos_nodes
V(g)$label <- ifelse(V(g)$name %in% union(top_pos_nodes, top_neg_nodes), V(g)$name, NA)
V(g)$label_sign <- case_when(
  V(g)$name %in% top_pos_nodes & V(g)$name %in% top_neg_nodes ~ "both",
  V(g)$name %in% top_pos_nodes ~ "positive",
  V(g)$name %in% top_neg_nodes ~ "negative",
  TRUE ~ NA_character_
)


# Filter weak edges
g_filt <- igraph::delete_edges(
  g,
  E(g)[abs(weight) < 0.05]
)
V(g_filt)$degree <- igraph::degree(g_filt)

V(g_filt)$size <- scales::rescale(V(g_filt)$degree, to = c(1, 12))

set.seed(42)

xy <- layout_components(g_filt, layout = layout_with_stress)
lay <- create_layout(g_filt, layout = xy)
lay_df <- as.data.frame(lay) %>% filter(!is.na(label))

V(g_filt)$deg <- igraph::degree(g_filt, mode = "all", loops = FALSE)
node_size <- scales::rescale(V(g_filt)$deg, to = c(3, 12))   # <--- AQUÍ se define

p <- ggraph(lay) +
  geom_edge_link(aes(color = sign,
                     edge_width = abs(weight),
                     edge_alpha = abs(weight)),
                 show.legend = TRUE) +
  scale_edge_color_manual(values = c(positive = "#1f78b4", negative = "#e31a1c"),
                          name = "Asociación") +
  scale_edge_width(range = c(0.2, 1.8)) +     
  scale_edge_alpha(range = c(0.5, 0.8)) +    
  guides(edge_width = "none", edge_alpha = "none") +
  
  geom_node_point(aes(fill = type),
                  shape = 21, size = node_size, color = "grey25", stroke = 0.4) +
  scale_fill_manual(values = c(Bacteria = "#E69F00", Fungi = "#0072B2", Unknown = "grey70"),
                    name = "Dominio") +
  theme_void() +
  coord_equal()

p #900*900


if (is.null(E(g)$sign)) {
  if (is.null(E(g)$weight)) stop("La red no tiene 'sign' ni 'weight' en las aristas.")
  E(g)$sign <- ifelse(E(g)$weight >= 0, "positive", "negative")
}

if (is.null(V(g)$type)) {
  V(g)$type <- ifelse(grepl("^bac_", V(g)$name), "Bacteria",
                      ifelse(grepl("^fun_", V(g)$name), "Fungi", "Unknown"))
}

edge_df <- igraph::as_data_frame(g, what = "edges") %>%
  mutate(
    from_type = V(g)$type[match(from, V(g)$name)],
    to_type   = V(g)$type[match(to,   V(g)$name)],
    # interaction type
    interaction = dplyr::case_when(
      from_type == "Bacteria" & to_type == "Bacteria" ~ "B-B",
      from_type == "Fungi"    & to_type == "Fungi"    ~ "F-F",
      TRUE                                            ~ "B-F"
    ),
    association = ifelse(as.character(sign) %in% c("positive","Positive"), "Positive", "Negative")
  )

prop_df <- edge_df %>%
  group_by(interaction, association) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(interaction) %>%
  mutate(proportion = count / sum(count))

edge_counts <- edge_df %>%
  count(interaction, association, name = "n") %>%
  group_by(interaction) %>%
  mutate(
    percentage = n / sum(n) * 100,
    label = paste0(n, "\n(", round(percentage, 1), "%)")
  ) %>%
  ungroup()

total_counts <- edge_counts %>%
  group_by(interaction) %>%
  summarise(total = sum(n), .groups = "drop")

edge_counts <- left_join(edge_counts, total_counts, by = "interaction")
assoc_colors <- c("Positive" = "#56B4E9",  
                  "Negative" = "#D55E00")  

#plot associations
y_max <- max(edge_counts$total, na.rm = TRUE)
neg_offset <- 0.05 * y_max  # 

ggplot(edge_counts, aes(x = interaction, y = n, fill = association)) +
  geom_bar(stat = "identity", color = "black") +
  
  
  geom_text(
    data = dplyr::filter(edge_counts, association == "Positive"),
    aes(label = label),
    position = position_stack(vjust = 0.5),
    color = "black", size = 4
  ) +
  
  geom_text(
    data = dplyr::filter(edge_counts, association == "Negative"),
    aes(y = total + neg_offset, label = label),
    color = "black", size = 4
  ) +
  
  scale_fill_manual(values = assoc_colors) +
  labs(x = "", y = "Number of Edges", fill = "Association") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    axis.line        = element_line(color = "black")
  )
#493*536

ggplot(edge_counts, aes(x = interaction, y = n, fill = association)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = assoc_colors) +
  labs(x = "", y = "Number of Edges", fill = "Association") +
  scale_y_continuous(limits = c(0, 200)) +   # <-- this controls Y-axis limit
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    axis.line        = element_line(color = "black")
  )
#493*536

#Obtain negative edges
negative_ff <- edge_df %>%
  filter(interaction == "F-F", association == "Negative")
negative_ff

negative_bf <- edge_df %>%
  filter(interaction == "B-F", association == "Negative")
negative_bf

positive_bf <- edge_df %>%
  filter(interaction == "B-F", association == "Positive")
positive_bf

negative_bb <- edge_df %>%
  filter(interaction == "B-B", association == "Negative")
negative_bb


}
  
net_BS_HC <- netConstruct(
  otu_BS_HC4,
  measure     = "spieceasi",
  normMethod  = "none",
  zeroMethod  = "pseudo",   
  filtTax     = "highestVar",
  filtTaxPar  = list(highestVar = fVar),
  
  # SpiecEasi parameters
  measurePar  = list(
    method           = "mb",        
    sel.criterion    = "stars",   
    pulsar.params    = list(
      rep.num          = 30,     
      thresh           = 0.05,     
      subsample.ratio  = 0.8,
      seed             = 2015
    ),
    nlambda          = 20,         
    lambda.min.ratio = 0.05        
  ),
  
  sparsMethod = "none",
  verbose     = 3,
  seed        = 2015
)

net_BS_HC_anal <- netAnalyze(net_BS_HC)
summary(net_BS_HC_anal)

edgelist <- net_BS_HC$edgelist1[
  order(net_BS_HC$edgelist1$adja, decreasing = TRUE), ]
edgelist |> head()
edgelist |> tail()

netcomi_graph <- SpiecEasi::adj2igraph(abs(net_BS_HC$adjaMat1))
netcomi_graph

str(net_BS_HC_anal)

edg <- net_BS_HC$edgelist1 %>%
  mutate(sign = ifelse(asso >= 0, "positive", "negative"),
         w_abs = abs(asso))

g <- graph_from_data_frame(edg, directed = FALSE)

V(g)$type <- ifelse(grepl("^bac_", V(g)$name), "Bacteria",
                    ifelse(grepl("^fun_", V(g)$name), "Fungi", "Unknown"))
V(g)$deg  <- igraph::degree(g) #degree
V(g)$eig  <- igraph::eigen_centrality(g)$vector #eigenvector centrality
V(g)$comm <- cluster_fast_greedy(g)$membership   # or cluster_louvain(g)

# Edge aesthetics
E(g)$weight <- edg$w_abs
E(g)$sign   <- edg$sign
E(g)$w_abs <- ifelse(is.na(E(g)$weight), 0, abs(E(g)$weight))

hub_n  <- 10
hub_ids <- order(V(g)$eig, decreasing = TRUE)[seq_len(min(hub_n, gorder(g)))]
V(g)$hub <- ifelse(seq_along(V(g)) %in% hub_ids, V(g)$name, NA)

if (is.null(V(g)$eig)) V(g)$eig <- igraph::evcent(g, scale = TRUE)$vector
if (is.null(V(g)$type)) {
  V(g)$type <- ifelse(grepl("^bac_", V(g)$name), "Bacteria",
                      ifelse(grepl("^fun_", V(g)$name), "Fungi", "Unknown"))
}
if (is.null(E(g)$sign)) {
  if (!is.null(E(g)$weight)) {
    E(g)$sign <- ifelse(E(g)$weight >= 0, "positive", "negative")
  } else stop("Falta 'sign' o 'weight' en edges para clasificar asociaciones.")
}

table(V(g)$type)


# ---------- Top 10 associations ----------
ed_df <- igraph::as_data_frame(g, what = "edges") %>%
  transmute(from, to, sign = as.character(sign))

count_by_sign <- ed_df %>%
  pivot_longer(c(from, to), names_to = "end", values_to = "node") %>%
  count(node, sign) %>%
  pivot_wider(names_from = sign, values_from = n, values_fill = 0) %>%
  rename(neg = negative, pos = positive)

top_pos_nodes <- count_by_sign %>% arrange(desc(pos)) %>% slice_head(n = 5) %>% pull(node)
top_neg_nodes <- count_by_sign %>% arrange(desc(neg)) %>% slice_head(n = 5) %>% pull(node)
top_pos_nodes
V(g)$label <- ifelse(V(g)$name %in% union(top_pos_nodes, top_neg_nodes), V(g)$name, NA)
V(g)$label_sign <- case_when(
  V(g)$name %in% top_pos_nodes & V(g)$name %in% top_neg_nodes ~ "both",
  V(g)$name %in% top_pos_nodes ~ "positive",
  V(g)$name %in% top_neg_nodes ~ "negative",
  TRUE ~ NA_character_
)


# Filter weak edges
g_filt <- igraph::delete_edges(
  g,
  E(g)[abs(weight) < 0.05]
)
V(g_filt)$degree <- igraph::degree(g_filt)

V(g_filt)$size <- scales::rescale(V(g_filt)$degree, to = c(1, 12))

set.seed(42)

xy <- layout_components(g_filt, layout = layout_with_stress)
lay <- create_layout(g_filt, layout = xy)
lay_df <- as.data.frame(lay) %>% filter(!is.na(label))

V(g_filt)$deg <- igraph::degree(g_filt, mode = "all", loops = FALSE)
node_size <- scales::rescale(V(g_filt)$deg, to = c(3, 12))   # <--- AQUÍ se define

p <- ggraph(lay) +
  geom_edge_link(aes(color = sign,
                     edge_width = abs(weight),
                     edge_alpha = abs(weight)),
                 show.legend = TRUE) +
  scale_edge_color_manual(values = c(positive = "#1f78b4", negative = "#e31a1c"),
                          name = "Asociación") +
  scale_edge_width(range = c(0.2, 1.8)) +     
  scale_edge_alpha(range = c(0.5, 0.8)) +    
  guides(edge_width = "none", edge_alpha = "none") +
  
  geom_node_point(aes(fill = type),
                  shape = 21, size = node_size, color = "grey25", stroke = 0.4) +
  scale_fill_manual(values = c(Bacteria = "#E69F00", Fungi = "#0072B2", Unknown = "grey70"),
                    name = "Dominio") +
  theme_void() +
  coord_equal()

p #900*900


if (is.null(E(g)$sign)) {
  if (is.null(E(g)$weight)) stop("La red no tiene 'sign' ni 'weight' en las aristas.")
  E(g)$sign <- ifelse(E(g)$weight >= 0, "positive", "negative")
}

if (is.null(V(g)$type)) {
  V(g)$type <- ifelse(grepl("^bac_", V(g)$name), "Bacteria",
                      ifelse(grepl("^fun_", V(g)$name), "Fungi", "Unknown"))
}

edge_df <- igraph::as_data_frame(g, what = "edges") %>%
  mutate(
    from_type = V(g)$type[match(from, V(g)$name)],
    to_type   = V(g)$type[match(to,   V(g)$name)],
    # interaction type
    interaction = dplyr::case_when(
      from_type == "Bacteria" & to_type == "Bacteria" ~ "B-B",
      from_type == "Fungi"    & to_type == "Fungi"    ~ "F-F",
      TRUE                                            ~ "B-F"
    ),
    association = ifelse(as.character(sign) %in% c("positive","Positive"), "Positive", "Negative")
  )

prop_df <- edge_df %>%
  group_by(interaction, association) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(interaction) %>%
  mutate(proportion = count / sum(count))

edge_counts <- edge_df %>%
  count(interaction, association, name = "n") %>%
  group_by(interaction) %>%
  mutate(
    percentage = n / sum(n) * 100,
    label = paste0(n, "\n(", round(percentage, 1), "%)")
  ) %>%
  ungroup()

total_counts <- edge_counts %>%
  group_by(interaction) %>%
  summarise(total = sum(n), .groups = "drop")

edge_counts <- left_join(edge_counts, total_counts, by = "interaction")
assoc_colors <- c("Positive" = "#56B4E9",  
                  "Negative" = "#D55E00")  

#plot associations
y_max <- max(edge_counts$total, na.rm = TRUE)
neg_offset <- 0.05 * y_max  # 

ggplot(edge_counts, aes(x = interaction, y = n, fill = association)) +
  geom_bar(stat = "identity", color = "black") +
  

  geom_text(
    data = dplyr::filter(edge_counts, association == "Positive"),
    aes(label = label),
    position = position_stack(vjust = 0.5),
    color = "black", size = 4
  ) +

  geom_text(
    data = dplyr::filter(edge_counts, association == "Negative"),
    aes(y = total + neg_offset, label = label),
    color = "black", size = 4
  ) +
  
  scale_fill_manual(values = assoc_colors) +
  labs(x = "", y = "Number of Edges", fill = "Association") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    axis.line        = element_line(color = "black")
  )
#493*536

ggplot(edge_counts, aes(x = interaction, y = n, fill = association)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = assoc_colors) +
  labs(x = "", y = "Number of Edges", fill = "Association") +
  scale_y_continuous(limits = c(0, 200)) +   # <-- this controls Y-axis limit
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    axis.line        = element_line(color = "black")
  )
#493*536

#Obtain negative edges
negative_ff <- edge_df %>%
  filter(interaction == "F-F", association == "Negative")
negative_ff

negative_bf <- edge_df %>%
  filter(interaction == "B-F", association == "Negative")
negative_bf

positive_bf <- edge_df %>%
  filter(interaction == "B-F", association == "Positive")
positive_bf

negative_bb <- edge_df %>%
  filter(interaction == "B-B", association == "Negative")
negative_bb


# Now for Soil surrounding roots
{
  fVar <- 250 
  net_SSR_UC <- netConstruct(
    otu_SSR_UC4,
    measure     = "spieceasi",
    normMethod  = "none",
    zeroMethod  = "pseudo",    
    filtTax     = "highestVar",
    filtTaxPar  = list(highestVar = fVar),
    
    # SpiecEasi parameters
    measurePar  = list(
      method           = "mb",       
      sel.criterion    = "stars",    
      pulsar.params    = list(
        rep.num          = 30,       
        thresh           = 0.05,      
        subsample.ratio  = 0.8,
        seed             = 2015
      ),
      nlambda          = 20,         
      lambda.min.ratio = 0.05         
    ),
    sparsMethod = "none",
    verbose     = 3,
    seed        = 2015
  )
  
  net_SSR_UC_anal <- netAnalyze(net_SSR_UC)
  summary(net_SSR_UC_anal)
  
  edgelist <- net_SSR_UC$edgelist1[
    order(net_SSR_UC$edgelist1$adja, decreasing = TRUE), ]
  edgelist |> head()
  edgelist |> tail()
  
  netcomi_graph <- SpiecEasi::adj2igraph(abs(net_SSR_UC$adjaMat1))
  netcomi_graph
  
  str(net_SSR_UC_anal)
  
  edg <- net_SSR_UC$edgelist1 %>%
    mutate(sign = ifelse(asso >= 0, "positive", "negative"),
           w_abs = abs(asso))
  
  g <- graph_from_data_frame(edg, directed = FALSE)
  
  V(g)$type <- ifelse(grepl("^bac_", V(g)$name), "Bacteria",
                      ifelse(grepl("^fun_", V(g)$name), "Fungi", "Unknown"))
  V(g)$deg  <- igraph::degree(g) #degree
  V(g)$eig  <- igraph::eigen_centrality(g)$vector #eigenvector centrality
  V(g)$comm <- cluster_fast_greedy(g)$membership   # or cluster_louvain(g)
  
  # Edge aesthetics
  E(g)$weight <- edg$w_abs
  E(g)$sign   <- edg$sign
  E(g)$w_abs <- ifelse(is.na(E(g)$weight), 0, abs(E(g)$weight))
  
  hub_n  <- 10
  hub_ids <- order(V(g)$eig, decreasing = TRUE)[seq_len(min(hub_n, gorder(g)))]
  V(g)$hub <- ifelse(seq_along(V(g)) %in% hub_ids, V(g)$name, NA)
  
  if (is.null(V(g)$eig)) V(g)$eig <- igraph::evcent(g, scale = TRUE)$vector
  if (is.null(V(g)$type)) {
    V(g)$type <- ifelse(grepl("^bac_", V(g)$name), "Bacteria",
                        ifelse(grepl("^fun_", V(g)$name), "Fungi", "Unknown"))
  }
  if (is.null(E(g)$sign)) {
    if (!is.null(E(g)$weight)) {
      E(g)$sign <- ifelse(E(g)$weight >= 0, "positive", "negative")
    } else stop("Falta 'sign' o 'weight' en edges para clasificar asociaciones.")
  }
  
  table(V(g)$type)
  
  
  # ---------- Top 10 associations ----------
  ed_df <- igraph::as_data_frame(g, what = "edges") %>%
    transmute(from, to, sign = as.character(sign))
  
  count_by_sign <- ed_df %>%
    pivot_longer(c(from, to), names_to = "end", values_to = "node") %>%
    count(node, sign) %>%
    pivot_wider(names_from = sign, values_from = n, values_fill = 0) %>%
    rename(neg = negative, pos = positive)
  
  top_pos_nodes <- count_by_sign %>% arrange(desc(pos)) %>% slice_head(n = 5) %>% pull(node)
  top_neg_nodes <- count_by_sign %>% arrange(desc(neg)) %>% slice_head(n = 5) %>% pull(node)
  top_pos_nodes
  V(g)$label <- ifelse(V(g)$name %in% union(top_pos_nodes, top_neg_nodes), V(g)$name, NA)
  V(g)$label_sign <- case_when(
    V(g)$name %in% top_pos_nodes & V(g)$name %in% top_neg_nodes ~ "both",
    V(g)$name %in% top_pos_nodes ~ "positive",
    V(g)$name %in% top_neg_nodes ~ "negative",
    TRUE ~ NA_character_
  )
  
  
  # Filter weak edges
  g_filt <- igraph::delete_edges(
    g,
    E(g)[abs(weight) < 0.05]
  )
  V(g_filt)$degree <- igraph::degree(g_filt)
  
  V(g_filt)$size <- scales::rescale(V(g_filt)$degree, to = c(1, 12))
  
  set.seed(42)
  
  xy <- layout_components(g_filt, layout = layout_with_stress)
  lay <- create_layout(g_filt, layout = xy)
  lay_df <- as.data.frame(lay) %>% filter(!is.na(label))
  
  V(g_filt)$deg <- igraph::degree(g_filt, mode = "all", loops = FALSE)
  node_size <- scales::rescale(V(g_filt)$deg, to = c(3, 12))   # <--- AQUÍ se define
  
  p <- ggraph(lay) +
    geom_edge_link(aes(color = sign,
                       edge_width = abs(weight),
                       edge_alpha = abs(weight)),
                   show.legend = TRUE) +
    scale_edge_color_manual(values = c(positive = "#1f78b4", negative = "#e31a1c"),
                            name = "Asociación") +
    scale_edge_width(range = c(0.2, 1.8)) +     
    scale_edge_alpha(range = c(0.5, 0.8)) +    
    guides(edge_width = "none", edge_alpha = "none") +
    
    geom_node_point(aes(fill = type),
                    shape = 21, size = node_size, color = "grey25", stroke = 0.4) +
    scale_fill_manual(values = c(Bacteria = "#E69F00", Fungi = "#0072B2", Unknown = "grey70"),
                      name = "Dominio") +
    theme_void() +
    coord_equal()
  
  p #900*900
  
  
  if (is.null(E(g)$sign)) {
    if (is.null(E(g)$weight)) stop("La red no tiene 'sign' ni 'weight' en las aristas.")
    E(g)$sign <- ifelse(E(g)$weight >= 0, "positive", "negative")
  }
  
  if (is.null(V(g)$type)) {
    V(g)$type <- ifelse(grepl("^bac_", V(g)$name), "Bacteria",
                        ifelse(grepl("^fun_", V(g)$name), "Fungi", "Unknown"))
  }
  
  edge_df <- igraph::as_data_frame(g, what = "edges") %>%
    mutate(
      from_type = V(g)$type[match(from, V(g)$name)],
      to_type   = V(g)$type[match(to,   V(g)$name)],
      # interaction type
      interaction = dplyr::case_when(
        from_type == "Bacteria" & to_type == "Bacteria" ~ "B-B",
        from_type == "Fungi"    & to_type == "Fungi"    ~ "F-F",
        TRUE                                            ~ "B-F"
      ),
      association = ifelse(as.character(sign) %in% c("positive","Positive"), "Positive", "Negative")
    )
  
  prop_df <- edge_df %>%
    group_by(interaction, association) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(interaction) %>%
    mutate(proportion = count / sum(count))
  
  edge_counts <- edge_df %>%
    count(interaction, association, name = "n") %>%
    group_by(interaction) %>%
    mutate(
      percentage = n / sum(n) * 100,
      label = paste0(n, "\n(", round(percentage, 1), "%)")
    ) %>%
    ungroup()
  
  total_counts <- edge_counts %>%
    group_by(interaction) %>%
    summarise(total = sum(n), .groups = "drop")
  
  edge_counts <- left_join(edge_counts, total_counts, by = "interaction")
  assoc_colors <- c("Positive" = "#56B4E9",  
                    "Negative" = "#D55E00")  
  
  #plot associations
  y_max <- max(edge_counts$total, na.rm = TRUE)
  neg_offset <- 0.05 * y_max  # 
  
  ggplot(edge_counts, aes(x = interaction, y = n, fill = association)) +
    geom_bar(stat = "identity", color = "black") +
    
    
    geom_text(
      data = dplyr::filter(edge_counts, association == "Positive"),
      aes(label = label),
      position = position_stack(vjust = 0.5),
      color = "black", size = 4
    ) +
    
    geom_text(
      data = dplyr::filter(edge_counts, association == "Negative"),
      aes(y = total + neg_offset, label = label),
      color = "black", size = 4
    ) +
    
    scale_fill_manual(values = assoc_colors) +
    labs(x = "", y = "Number of Edges", fill = "Association") +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA),
      axis.line        = element_line(color = "black")
    )
  #493*536
  
  ggplot(edge_counts, aes(x = interaction, y = n, fill = association)) +
    geom_bar(stat = "identity", color = "black") +
    scale_fill_manual(values = assoc_colors) +
    labs(x = "", y = "Number of Edges", fill = "Association") +
    scale_y_continuous(limits = c(0, 200)) +   # <-- this controls Y-axis limit
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA),
      axis.line        = element_line(color = "black")
    )
  #493*536
  
  #Obtain negative edges
  negative_ff <- edge_df %>%
    filter(interaction == "F-F", association == "Negative")
  negative_ff
  
  negative_bf <- edge_df %>%
    filter(interaction == "B-F", association == "Negative")
  negative_bf
  
  positive_bf <- edge_df %>%
    filter(interaction == "B-F", association == "Positive")
  positive_bf
  
  negative_bb <- edge_df %>%
    filter(interaction == "B-B", association == "Negative")
  negative_bb
  
}


#BS MC

net_SSR_MC <- netConstruct(
  otu_SSR_MC4,
  measure     = "spieceasi",
  normMethod  = "none",
  zeroMethod  = "pseudo",   
  filtTax     = "highestVar",
  filtTaxPar  = list(highestVar = fVar),
  
  # SpiecEasi parameters
  measurePar  = list(
    method           = "mb",      
    sel.criterion    = "stars",   
    pulsar.params    = list(
      rep.num          = 30,       
      thresh           = 0.05,      
      subsample.ratio  = 0.8,
      seed             = 2015
    ),
    nlambda          = 20,         
    lambda.min.ratio = 0.05        
  ),
  sparsMethod = "none",
  verbose     = 3,
  seed        = 2015
)

net_SSR_MC_anal <- netAnalyze(net_SSR_MC)
summary(net_SSR_MC_anal)

edgelist <- net_SSR_MC$edgelist1[
  order(net_SSR_MC$edgelist1$adja, decreasing = TRUE), ]
edgelist |> head()
edgelist |> tail()

netcomi_graph <- SpiecEasi::adj2igraph(abs(net_SSR_MC$adjaMat1))
netcomi_graph

str(net_SSR_MC_anal)

edg <- net_SSR_MC$edgelist1 %>%
  mutate(sign = ifelse(asso >= 0, "positive", "negative"),
         w_abs = abs(asso))

g <- graph_from_data_frame(edg, directed = FALSE)

V(g)$type <- ifelse(grepl("^bac_", V(g)$name), "Bacteria",
                    ifelse(grepl("^fun_", V(g)$name), "Fungi", "Unknown"))
V(g)$deg  <- igraph::degree(g) #degree
V(g)$eig  <- igraph::eigen_centrality(g)$vector #eigenvector centrality
V(g)$comm <- cluster_fast_greedy(g)$membership   # or cluster_louvain(g)

# Edge aesthetics
E(g)$weight <- edg$w_abs
E(g)$sign   <- edg$sign
E(g)$w_abs <- ifelse(is.na(E(g)$weight), 0, abs(E(g)$weight))

hub_n  <- 10
hub_ids <- order(V(g)$eig, decreasing = TRUE)[seq_len(min(hub_n, gorder(g)))]
V(g)$hub <- ifelse(seq_along(V(g)) %in% hub_ids, V(g)$name, NA)

if (is.null(V(g)$eig)) V(g)$eig <- igraph::evcent(g, scale = TRUE)$vector
if (is.null(V(g)$type)) {
  V(g)$type <- ifelse(grepl("^bac_", V(g)$name), "Bacteria",
                      ifelse(grepl("^fun_", V(g)$name), "Fungi", "Unknown"))
}
if (is.null(E(g)$sign)) {
  if (!is.null(E(g)$weight)) {
    E(g)$sign <- ifelse(E(g)$weight >= 0, "positive", "negative")
  } else stop("Falta 'sign' o 'weight' en edges para clasificar asociaciones.")
}

table(V(g)$type)


# ---------- Top 10 associations ----------
ed_df <- igraph::as_data_frame(g, what = "edges") %>%
  transmute(from, to, sign = as.character(sign))

count_by_sign <- ed_df %>%
  pivot_longer(c(from, to), names_to = "end", values_to = "node") %>%
  count(node, sign) %>%
  pivot_wider(names_from = sign, values_from = n, values_fill = 0) %>%
  rename(neg = negative, pos = positive)

top_pos_nodes <- count_by_sign %>% arrange(desc(pos)) %>% slice_head(n = 5) %>% pull(node)
top_neg_nodes <- count_by_sign %>% arrange(desc(neg)) %>% slice_head(n = 5) %>% pull(node)
top_pos_nodes
V(g)$label <- ifelse(V(g)$name %in% union(top_pos_nodes, top_neg_nodes), V(g)$name, NA)
V(g)$label_sign <- case_when(
  V(g)$name %in% top_pos_nodes & V(g)$name %in% top_neg_nodes ~ "both",
  V(g)$name %in% top_pos_nodes ~ "positive",
  V(g)$name %in% top_neg_nodes ~ "negative",
  TRUE ~ NA_character_
)


# Filter weak edges
g_filt <- igraph::delete_edges(
  g,
  E(g)[abs(weight) < 0.05]
)
V(g_filt)$degree <- igraph::degree(g_filt)

V(g_filt)$size <- scales::rescale(V(g_filt)$degree, to = c(1, 12))

set.seed(42)

xy <- layout_components(g_filt, layout = layout_with_stress)
lay <- create_layout(g_filt, layout = xy)
lay_df <- as.data.frame(lay) %>% filter(!is.na(label))

V(g_filt)$deg <- igraph::degree(g_filt, mode = "all", loops = FALSE)
node_size <- scales::rescale(V(g_filt)$deg, to = c(3, 12))   # <--- AQUÍ se define

p <- ggraph(lay) +
  geom_edge_link(aes(color = sign,
                     edge_width = abs(weight),
                     edge_alpha = abs(weight)),
                 show.legend = TRUE) +
  scale_edge_color_manual(values = c(positive = "#1f78b4", negative = "#e31a1c"),
                          name = "Asociación") +
  scale_edge_width(range = c(0.2, 1.8)) +     
  scale_edge_alpha(range = c(0.5, 0.8)) +    
  guides(edge_width = "none", edge_alpha = "none") +
  
  geom_node_point(aes(fill = type),
                  shape = 21, size = node_size, color = "grey25", stroke = 0.4) +
  scale_fill_manual(values = c(Bacteria = "#E69F00", Fungi = "#0072B2", Unknown = "grey70"),
                    name = "Dominio") +
  theme_void() +
  coord_equal()

p #900*900


if (is.null(E(g)$sign)) {
  if (is.null(E(g)$weight)) stop("La red no tiene 'sign' ni 'weight' en las aristas.")
  E(g)$sign <- ifelse(E(g)$weight >= 0, "positive", "negative")
}

if (is.null(V(g)$type)) {
  V(g)$type <- ifelse(grepl("^bac_", V(g)$name), "Bacteria",
                      ifelse(grepl("^fun_", V(g)$name), "Fungi", "Unknown"))
}

edge_df <- igraph::as_data_frame(g, what = "edges") %>%
  mutate(
    from_type = V(g)$type[match(from, V(g)$name)],
    to_type   = V(g)$type[match(to,   V(g)$name)],
    # interaction type
    interaction = dplyr::case_when(
      from_type == "Bacteria" & to_type == "Bacteria" ~ "B-B",
      from_type == "Fungi"    & to_type == "Fungi"    ~ "F-F",
      TRUE                                            ~ "B-F"
    ),
    association = ifelse(as.character(sign) %in% c("positive","Positive"), "Positive", "Negative")
  )

prop_df <- edge_df %>%
  group_by(interaction, association) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(interaction) %>%
  mutate(proportion = count / sum(count))

edge_counts <- edge_df %>%
  count(interaction, association, name = "n") %>%
  group_by(interaction) %>%
  mutate(
    percentage = n / sum(n) * 100,
    label = paste0(n, "\n(", round(percentage, 1), "%)")
  ) %>%
  ungroup()

total_counts <- edge_counts %>%
  group_by(interaction) %>%
  summarise(total = sum(n), .groups = "drop")

edge_counts <- left_join(edge_counts, total_counts, by = "interaction")
assoc_colors <- c("Positive" = "#56B4E9",  
                  "Negative" = "#D55E00")  

#plot associations
y_max <- max(edge_counts$total, na.rm = TRUE)
neg_offset <- 0.05 * y_max  # 

ggplot(edge_counts, aes(x = interaction, y = n, fill = association)) +
  geom_bar(stat = "identity", color = "black") +
  
  
  geom_text(
    data = dplyr::filter(edge_counts, association == "Positive"),
    aes(label = label),
    position = position_stack(vjust = 0.5),
    color = "black", size = 4
  ) +
  
  geom_text(
    data = dplyr::filter(edge_counts, association == "Negative"),
    aes(y = total + neg_offset, label = label),
    color = "black", size = 4
  ) +
  
  scale_fill_manual(values = assoc_colors) +
  labs(x = "", y = "Number of Edges", fill = "Association") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    axis.line        = element_line(color = "black")
  )
#493*536

ggplot(edge_counts, aes(x = interaction, y = n, fill = association)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = assoc_colors) +
  labs(x = "", y = "Number of Edges", fill = "Association") +
  scale_y_continuous(limits = c(0, 200)) +   # <-- this controls Y-axis limit
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    axis.line        = element_line(color = "black")
  )
#493*536

#Obtain negative edges
negative_ff <- edge_df %>%
  filter(interaction == "F-F", association == "Negative")
negative_ff

negative_bf <- edge_df %>%
  filter(interaction == "B-F", association == "Negative")
negative_bf

positive_bf <- edge_df %>%
  filter(interaction == "B-F", association == "Positive")
positive_bf

negative_bb <- edge_df %>%
  filter(interaction == "B-B", association == "Negative")
negative_bb




#SSR_HC
net_SSR_HC <- netConstruct(
  otu_SSR_HC4,
  measure     = "spieceasi",
  normMethod  = "none",
  zeroMethod  = "pseudo",    
  filtTax     = "highestVar",
  filtTaxPar  = list(highestVar = fVar),
  
  # Parámetros de SpiecEasi
  measurePar  = list(
    method           = "mb",        
    sel.criterion    = "stars",    
    pulsar.params    = list(
      rep.num          = 30,       
      thresh           = 0.05,      
      subsample.ratio  = 0.8,
      seed             = 2015
    ),
    nlambda          = 20,          
    lambda.min.ratio = 0.05         
  ),
  
  sparsMethod = "none",
  verbose     = 3,
  seed        = 2015
)

net_SSR_HC_anal <- netAnalyze(net_SSR_HC)
summary(net_SSR_HC_anal)

edgelist <- net_SSR_HC$edgelist1[
  order(net_SSR_HC$edgelist1$adja, decreasing = TRUE), ]
edgelist |> head()
edgelist |> tail()

netcomi_graph <- SpiecEasi::adj2igraph(abs(net_SSR_HC$adjaMat1))
netcomi_graph

str(net_SSR_HC_anal)

edg <- net_SSR_HC$edgelist1 %>%
  mutate(sign = ifelse(asso >= 0, "positive", "negative"),
         w_abs = abs(asso))

g <- graph_from_data_frame(edg, directed = FALSE)

V(g)$type <- ifelse(grepl("^bac_", V(g)$name), "Bacteria",
                    ifelse(grepl("^fun_", V(g)$name), "Fungi", "Unknown"))
V(g)$deg  <- igraph::degree(g) #degree
V(g)$eig  <- igraph::eigen_centrality(g)$vector #eigenvector centrality
V(g)$comm <- cluster_fast_greedy(g)$membership   # or cluster_louvain(g)

# Edge aesthetics
E(g)$weight <- edg$w_abs
E(g)$sign   <- edg$sign
E(g)$w_abs <- ifelse(is.na(E(g)$weight), 0, abs(E(g)$weight))

hub_n  <- 10
hub_ids <- order(V(g)$eig, decreasing = TRUE)[seq_len(min(hub_n, gorder(g)))]
V(g)$hub <- ifelse(seq_along(V(g)) %in% hub_ids, V(g)$name, NA)

if (is.null(V(g)$eig)) V(g)$eig <- igraph::evcent(g, scale = TRUE)$vector
if (is.null(V(g)$type)) {
  V(g)$type <- ifelse(grepl("^bac_", V(g)$name), "Bacteria",
                      ifelse(grepl("^fun_", V(g)$name), "Fungi", "Unknown"))
}
if (is.null(E(g)$sign)) {
  if (!is.null(E(g)$weight)) {
    E(g)$sign <- ifelse(E(g)$weight >= 0, "positive", "negative")
  } else stop("Falta 'sign' o 'weight' en edges para clasificar asociaciones.")
}

table(V(g)$type)


# ---------- Top 10 associations ----------
ed_df <- igraph::as_data_frame(g, what = "edges") %>%
  transmute(from, to, sign = as.character(sign))

count_by_sign <- ed_df %>%
  pivot_longer(c(from, to), names_to = "end", values_to = "node") %>%
  count(node, sign) %>%
  pivot_wider(names_from = sign, values_from = n, values_fill = 0) %>%
  rename(neg = negative, pos = positive)

top_pos_nodes <- count_by_sign %>% arrange(desc(pos)) %>% slice_head(n = 5) %>% pull(node)
top_neg_nodes <- count_by_sign %>% arrange(desc(neg)) %>% slice_head(n = 5) %>% pull(node)
top_pos_nodes
V(g)$label <- ifelse(V(g)$name %in% union(top_pos_nodes, top_neg_nodes), V(g)$name, NA)
V(g)$label_sign <- case_when(
  V(g)$name %in% top_pos_nodes & V(g)$name %in% top_neg_nodes ~ "both",
  V(g)$name %in% top_pos_nodes ~ "positive",
  V(g)$name %in% top_neg_nodes ~ "negative",
  TRUE ~ NA_character_
)


# Filter weak edges
g_filt <- igraph::delete_edges(
  g,
  E(g)[abs(weight) < 0.05]
)
V(g_filt)$degree <- igraph::degree(g_filt)

V(g_filt)$size <- scales::rescale(V(g_filt)$degree, to = c(1, 12))

set.seed(42)

xy <- layout_components(g_filt, layout = layout_with_stress)
lay <- create_layout(g_filt, layout = xy)
lay_df <- as.data.frame(lay) %>% filter(!is.na(label))

V(g_filt)$deg <- igraph::degree(g_filt, mode = "all", loops = FALSE)
node_size <- scales::rescale(V(g_filt)$deg, to = c(3, 12))   # <--- AQUÍ se define

p <- ggraph(lay) +
  geom_edge_link(aes(color = sign,
                     edge_width = abs(weight),
                     edge_alpha = abs(weight)),
                 show.legend = TRUE) +
  scale_edge_color_manual(values = c(positive = "#1f78b4", negative = "#e31a1c"),
                          name = "Asociación") +
  scale_edge_width(range = c(0.2, 1.8)) +     
  scale_edge_alpha(range = c(0.5, 0.8)) +    
  guides(edge_width = "none", edge_alpha = "none") +
  
  geom_node_point(aes(fill = type),
                  shape = 21, size = node_size, color = "grey25", stroke = 0.4) +
  scale_fill_manual(values = c(Bacteria = "#E69F00", Fungi = "#0072B2", Unknown = "grey70"),
                    name = "Dominio") +
  theme_void() +
  coord_equal()

p #900*900


if (is.null(E(g)$sign)) {
  if (is.null(E(g)$weight)) stop("La red no tiene 'sign' ni 'weight' en las aristas.")
  E(g)$sign <- ifelse(E(g)$weight >= 0, "positive", "negative")
}

if (is.null(V(g)$type)) {
  V(g)$type <- ifelse(grepl("^bac_", V(g)$name), "Bacteria",
                      ifelse(grepl("^fun_", V(g)$name), "Fungi", "Unknown"))
}

edge_df <- igraph::as_data_frame(g, what = "edges") %>%
  mutate(
    from_type = V(g)$type[match(from, V(g)$name)],
    to_type   = V(g)$type[match(to,   V(g)$name)],
    # interaction type
    interaction = dplyr::case_when(
      from_type == "Bacteria" & to_type == "Bacteria" ~ "B-B",
      from_type == "Fungi"    & to_type == "Fungi"    ~ "F-F",
      TRUE                                            ~ "B-F"
    ),
    association = ifelse(as.character(sign) %in% c("positive","Positive"), "Positive", "Negative")
  )

prop_df <- edge_df %>%
  group_by(interaction, association) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(interaction) %>%
  mutate(proportion = count / sum(count))

edge_counts <- edge_df %>%
  count(interaction, association, name = "n") %>%
  group_by(interaction) %>%
  mutate(
    percentage = n / sum(n) * 100,
    label = paste0(n, "\n(", round(percentage, 1), "%)")
  ) %>%
  ungroup()

total_counts <- edge_counts %>%
  group_by(interaction) %>%
  summarise(total = sum(n), .groups = "drop")

edge_counts <- left_join(edge_counts, total_counts, by = "interaction")
assoc_colors <- c("Positive" = "#56B4E9",  
                  "Negative" = "#D55E00")  

#plot associations
y_max <- max(edge_counts$total, na.rm = TRUE)
neg_offset <- 0.05 * y_max  # 

ggplot(edge_counts, aes(x = interaction, y = n, fill = association)) +
  geom_bar(stat = "identity", color = "black") +
  
  
  geom_text(
    data = dplyr::filter(edge_counts, association == "Positive"),
    aes(label = label),
    position = position_stack(vjust = 0.5),
    color = "black", size = 4
  ) +
  
  geom_text(
    data = dplyr::filter(edge_counts, association == "Negative"),
    aes(y = total + neg_offset, label = label),
    color = "black", size = 4
  ) +
  
  scale_fill_manual(values = assoc_colors) +
  labs(x = "", y = "Number of Edges", fill = "Association") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    axis.line        = element_line(color = "black")
  )
#493*536

ggplot(edge_counts, aes(x = interaction, y = n, fill = association)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = assoc_colors) +
  labs(x = "", y = "Number of Edges", fill = "Association") +
  scale_y_continuous(limits = c(0, 200)) +   # <-- this controls Y-axis limit
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    axis.line        = element_line(color = "black")
  )
#493*536

#Obtain negative edges
negative_ff <- edge_df %>%
  filter(interaction == "F-F", association == "Negative")
negative_ff

negative_bf <- edge_df %>%
  filter(interaction == "B-F", association == "Negative")
negative_bf

positive_bf <- edge_df %>%
  filter(interaction == "B-F", association == "Positive")
positive_bf

negative_bb <- edge_df %>%
  filter(interaction == "B-B", association == "Negative")
negative_bb






#Summaries
summary(net_BS_UC_anal)
summary(net_BS_MC_anal)
summary(net_BS_HC_anal)
summary(net_SSR_UC_anal)
summary(net_SSR_MC_anal)
summary(net_SSR_HC_anal)

