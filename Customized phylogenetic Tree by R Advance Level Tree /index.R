# Load necessary libraries1
library(ape)
library(phangorn)
library(seqinr)
library(msa)
library(ggtree)
library(ggplot2)
library(phytools) # This library contains the midpoint.root function
library(ggimage)


if (!requireNamespace("treeio", quietly = TRUE)) {
  install.packages("treeio")
}
library(treeio)


BiocManager::install("msa", force = TRUE)
BiocManager::install("ggtree")


# Read your protein sequences
protein_seqs <- readAAStringSet("your_protein_sequences.fasta")


# Perform multiple sequence alignment
protein_msa <- msa(protein_seqs , method = "ClustalW")


# Convert MSA to phyDat
protein_phyDat <- as.phyDat(protein_msa, type = "AA", model = "LG")


# Calculate distance matrix (choose an appropriate model)
dist_matrix <- dist.ml(protein_phyDat)


# Construct a tree (example using Neighbor-Joining)
protein_tree <- nj(dist_matrix)


# Visualize the tree
plot(protein_tree)


# Optionally, reroot the tree (e.g., using the midpoint)
rerooted_tree <- midpoint.root(protein_tree)

##Print the tree to verify its content
cat("Tree in Newick format:\n")
cat(write.tree(rerooted_tree), "\n")


# Save the tree in Newick format
write.tree(rerooted_tree, file = "protein_tree.newick")

# Save the tree in Nexus format
write.nexus(rerooted_tree, file = "protein_tree.nexus")


 # Save the tree in Newick format using treeio
tree_data <- as.treedata(rerooted_tree)
write.tree(tree_data, file = "protein_tree_treeio.newick")


# Convert tree to ggtree object for advanced visualization
ggtree_tree <- as.phylo(rerooted_tree)

ggtree_tree
plot(ggtree_tree)


## advance code for styling tree
# Plot the tree with advanced styling
ggtree_plot <- ggtree(ggtree_tree) +
  geom_tiplab(aes(label=label), size=3, color="blue", hjust=-0.2) +
  geom_text(aes(label=node), vjust=-0.5, size=2, color="red") +
  geom_treescale(width=0.1, x=0, y=-2, color="darkgreen") +
  theme_tree2() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title="Phylogenetic Tree")

# Display the tree
print(ggtree_plot)

# Save the plot to a file
ggsave("protein_tree_only_labels.png", ggtree_plot, width=10, height=8)



### Circular Tree 

# Plot the tree with advanced styling in circular layout
ggtree_plot <- ggtree(ggtree_tree, layout="circular", size = 2) +
  geom_tiplab(aes(label=label), size=3, color="blue", align=TRUE, linesize=0.5) +
  geom_text2(aes(subset = !isTip, label = node), vjust=-0.5, size=3, color="red") +
  geom_treescale(width=0.1, color="darkgreen") +
  theme_tree2() +
  theme(
    plot.title = element_text(hjust = 0, size=16, face="bold", color = "pink", margin=margin(b=100)),
    plot.margin = margin(10, 100, 100, 100)  # Add margins around the plot
  ) +
  labs(title="Circular Phylogenetic Tree 😍 ")

# Display the tree
print(ggtree_plot)

# Save the plot to a file
ggsave("protein_tree_circularr.png", ggtree_plot, width=10, height=10)



##Fan Tree  😍

# Fan tree
fan_plot <- ggtree(ggtree_tree, layout="fan") +
  geom_tiplab(aes(label=label), size=4, color="blue", align=TRUE, linesize=0.5) +
  geom_text2(aes(subset = !isTip, label = node), vjust=-0.5, size=3, color="red") +
  geom_treescale(width=0.1, color="darkgreen") +
  theme_tree2() +
  theme(
    plot.title = element_text(hjust = 0.5, size=16, face="bold", margin=margin(b=20)),
    plot.margin = margin(20, 20, 20, 20)
  ) +
  labs(title="😍 Fan Phylogenetic Tree 😍")

# Display the tree
print(fan_plot)

# Save the plot to a file
ggsave("fan_tree.png", fan_plot, width=12, height=12)


## Unrooted Tree

# Unrooted tree
unrooted_plot <- ggtree(ggtree_tree, layout="unrooted") +
  geom_tiplab(aes(label=label), size=3, color="blue", hjust=-0.2) +
  geom_text2(aes(subset = !isTip, label = node), vjust=-0.5, size=2, color="red") +
  theme_tree2() +
  theme(
    plot.title = element_text(hjust = 0.5, size=16, face="bold", margin=margin(b=20)),
    plot.margin = margin(20, 20, 20, 20)
  ) +
  labs(title="😍 Unrooted Phylogenetic Tree 😍 ")

# Display the tree
print(unrooted_plot)

# Save the plot to a file
ggsave("unrooted_tree.png", unrooted_plot, width=12, height=12)


## Radial Tree

# Radial tree
radial_plot <- ggtree(ggtree_tree, layout="radial", size = 1) +
  geom_tiplab(aes(label=label), size=4, color="blue", align=TRUE, linesize=0.5) +
  geom_text2(aes(subset = !isTip, label = node), vjust=-0.5, size=2, color="red") +
  geom_treescale(width=0.1, color="darkgreen") +
  theme_tree2() +
  theme(
    plot.title = element_text(hjust = 0.5, size=16, face="bold", margin=margin(b=20) , color= "#f5007a"),
    plot.margin = margin(70, 70, 70, 70)
  ) +
  labs(title="😍 Radial Phylogenetic Tree 😍 ")

# Display the tree
print(radial_plot)

# Save the plot to a file
ggsave("radial_tree.png", radial_plot, width=12, height=12)






