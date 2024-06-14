library(ape)
library(phangorn)
library(seqinr)
library(msa)

BiocManager::install("msa")



# Read your protein sequences
protein_seqs <- readAAStringSet("your_protein_sequences.fasta")


# Perform multiple sequence alignment
protein_msa <- msa(protein_seqs)


# Convert MSA to phyDat
protein_phyDat <- as.phyDat(protein_msa, type = "AA", model = "WAG")


# Calculate distance matrix (choose an appropriate model)
dist_matrix <- dist.ml(protein_phyDat)


# Construct a tree (example using Neighbor-Joining)
protein_tree <- nj(dist_matrix)


# Visualize the tree
plot(protein_tree)
